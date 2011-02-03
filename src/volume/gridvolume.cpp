/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/volume.h>
#include <mitsuba/core/fstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

#if defined(__LINUX__) || defined(__OSX__)
#include <sys/mman.h>
#include <fcntl.h>
#endif

MTS_NAMESPACE_BEGIN

/**
 * This class implements a simple binary exchange format
 * for discretized volume data.
 */
class GridDataSource : public VolumeDataSource {
public:
	GridDataSource(const Properties &props) 
		: VolumeDataSource(props), m_overrideAABB(false) {
		std::string filename = props.getString("filename");
		loadFromFile(props.getString("filename"));
		m_volumeToWorld = props.getTransform("toWorld", Transform());

		if (props.hasProperty("min") && props.hasProperty("max")) {
			/* Optionally allow to use an AABB other than 
			   the one specified by the grid file */
			m_originalAABB.min = props.getPoint("min");
			m_originalAABB.max = props.getPoint("max");
			m_overrideAABB = true;
		}

		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_volumeToWorld(m_originalAABB.getCorner(i)));

		/**
		* When 'sendData' is set to false, only the filename 
		* is transmitted. A following unserialization of the 
		* stream causes the implementation to then look for 
		* the file (which had better exist if unserialization 
		* occurs on a remote machine)
		*/
		m_sendData = props.getBoolean("sendData", true);
	}

	GridDataSource(Stream *stream, InstanceManager *manager) 
	: VolumeDataSource(stream, manager) {
		m_volumeToWorld = Transform(stream);
		m_originalAABB = AABB(stream);
		m_sendData = stream->readBool();
		m_overrideAABB = stream->readBool();
		if (m_sendData) { 
			m_res = Vector3i(stream);
			m_channels = stream->readInt();
			m_filename = stream->readString();
			size_t nEntries = m_res.x*m_res.y*m_res.z*m_channels;
			m_data = new float[nEntries];
			stream->readSingleArray(m_data, nEntries);
			m_fromStream = true;
		} else {
			std::string filename = stream->readString();
			loadFromFile(filename);
		}
		configure();
	}

	virtual ~GridDataSource() {
#if defined(__LINUX__) || defined(__OSX__)
		if (m_fromStream) {
			delete[] m_data;
		} else {
			Log(EDebug, "Unmapping \"%s\" from memory", 
				m_filename.file_string().c_str()); 
			int retval = munmap(m_mmapPtr, m_mmapSize);
			if (retval != 0)
				Log(EError, "munmap(): unable to unmap memory!");
		}
#elif defined(WIN32)
		if (!UnmapViewOfFile(m_mmapPtr))
			Log(EError, "UnmapViewOfFile(): unable to unmap memory: %s", lastErrorText().c_str());
		if (!CloseHandle(m_fileMapping))
			Log(EError, "CloseHandle(): unable to close file mapping: %s", lastErrorText().c_str());
		if (!CloseHandle(m_file))
			Log(EError, "CloseHandle(): unable to close file: %s", lastErrorText().c_str());
#else
		delete[] m_data;
#endif
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);

		m_volumeToWorld.serialize(stream);
		m_originalAABB.serialize(stream);
		stream->writeBool(m_sendData);
		stream->writeBool(m_overrideAABB);

		if (m_sendData) {
			m_res.serialize(stream);
			stream->writeInt(m_channels);
			stream->writeString(m_filename.file_string());
			size_t nEntries = m_res.x*m_res.y*m_res.z*m_channels;
			stream->writeSingleArray(m_data, nEntries);
		} else {
			stream->writeString(m_filename.file_string());
		}
	}
		
	void configure() {
		Vector extents(m_originalAABB.getExtents());
		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(
				(m_res[0] - 1) / extents[0],
				(m_res[1] - 1) / extents[1],
				(m_res[2] - 1) / extents[2])
			) * Transform::translate(-Vector(m_originalAABB.min)) * m_worldToVolume;
		m_stepSize = std::numeric_limits<Float>::infinity();
		for (int i=0; i<3; ++i)
			m_stepSize = std::min(m_stepSize, extents[i] / (Float) (m_res[i]-1));
	}

	void loadFromFile(const std::string &filename) {
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(filename);
		ref<FileStream> stream = new FileStream(resolved, FileStream::EReadOnly);
		stream->setByteOrder(Stream::ELittleEndian);
		m_filename = filename;
		m_fromStream = false;

		char header[3];
		stream->read(header, 3);
		if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
			Log(EError, "Encountered an invalid volume data file (incorrect header identifier)");
		uint8_t version;
		stream->read(&version, 1);
		if (version != 3)
			Log(EError, "Encountered an invalid volume data file (incorrect file version)");
		int type = stream->readInt();
		if (type != 1)
			Log(EError, "Encountered an invalid volume data file (incorrect data type)");

		int xres = stream->readInt(), yres=stream->readInt(), zres=stream->readInt();
		m_res = Vector3i(xres, yres, zres);
		m_channels = stream->readInt();

		size_t nEntries = m_res.x*m_res.y*m_res.z*m_channels;
		Float xmin = stream->readSingle(), ymin = stream->readSingle(), zmin = stream->readSingle();
		Float xmax = stream->readSingle(), ymax = stream->readSingle(), zmax = stream->readSingle();
		if (!m_overrideAABB)
			m_originalAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));

#if defined(__LINUX__) || defined(__OSX__)
		Log(EDebug, "Mapping \"%s\" into memory: %ix%ix%i (%i channels), %i KiB, %s", filename.c_str(), 
			m_res.x, m_res.y, m_res.z, m_channels, nEntries*sizeof(float)/1024,
			m_originalAABB.toString().c_str());
		stream->close();
		int fd = open(resolved.file_string().c_str(), O_RDONLY);
		if (fd == -1)
			Log(EError, "Could not open \"%s\"!", m_filename.file_string().c_str());
		m_mmapSize = (nEntries+12)*sizeof(float);
		m_mmapPtr = mmap(NULL, m_mmapSize, PROT_READ, MAP_SHARED, fd, 0);
		if (m_mmapPtr == NULL)
			Log(EError, "Could not map \"%s\" to memory!", m_filename.file_string().c_str());
		m_data = ((float *) m_mmapPtr) + 12;
		if (close(fd) != 0)
			Log(EError, "close(): unable to close file!");
#elif defined(WIN32)
		Log(EDebug, "Mapping \"%s\" into memory: %ix%ix%i (%i channels), %i KiB, %s", filename.c_str(), 
			m_res.x, m_res.y, m_res.z, m_channels, nEntries*sizeof(float)/1024,
			m_originalAABB.toString().c_str());
		stream->close();
		m_file = CreateFile(resolved.file_string().c_str(), GENERIC_READ, 
			FILE_SHARE_READ, NULL, OPEN_EXISTING, 
			FILE_ATTRIBUTE_NORMAL, NULL);
		if (m_file == INVALID_HANDLE_VALUE)
			Log(EError, "Could not open \"%s\": %s", m_filename.file_string().c_str(),
				lastErrorText().c_str());
		m_fileMapping = CreateFileMapping(m_file, NULL, PAGE_READONLY, 0, 0, NULL);
		if (m_fileMapping == NULL)
			Log(EError, "CreateFileMapping: Could not map \"%s\" to memory: %s", 
				m_filename.file_string().c_str(), lastErrorText().c_str());
		m_mmapPtr = (float *) MapViewOfFile(m_fileMapping, FILE_MAP_READ, 0, 0, 0);
		if (m_mmapPtr == NULL)
			Log(EError, "MapViewOfFile: Could not map \"%s\" to memory: %s", 
				m_filename.file_string().c_str(), lastErrorText().c_str());
		m_data = ((float *) m_mmapPtr) + 12;
#else
		Log(EInfo, "Loading \"%s\": %ix%ix%i (%i channels), %i KiB, %s", filename.c_str(), 
			m_res.x, m_res.y, m_res.z, m_channels, nEntries*sizeof(float)/1024,
			m_originalAABB.toString().c_str());
		m_data = new float[nEntries];
		stream->read(m_data, nEntries*sizeof(float));
		stream->close();
#endif
	}

	/**
	 * This is needed since Mitsuba might be 
	 * compiled with either single/double precision
	 */
	struct float3 {
		float value[3];

		inline float3(float a, float b, float c) {
			value[0] = a; value[1] = b; value[2] = c;
		}

		inline float3(float *x) {
			value[0] = x[0]; value[1] = x[1]; value[1] = x[2];
		}

		inline float3 operator*(Float v) const {
			return float3(value[0]*v, value[1]*v, value[2]*v);
		}
		
		inline float3 operator+(const float3 &f2) const {
			return float3(value[0]+f2.value[0], value[1]+f2.value[1], value[2]+f2.value[2]);
		}

		inline Spectrum toSpectrum() const {
			Spectrum result;
			result.fromLinearRGB(value[0], value[1], value[2]);
			return result;
		}
		
		inline Vector toVector() const {
			return Vector(value[0], value[1], value[2]);
		}
	};

	Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		if (p.x < 0 || p.y < 0 || p.z < 0)
			return 0.0f;
		const int x1 = (int) p.x, y1 = (int) p.y, z1 = (int) p.z,
					x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x || 
			y2 >= m_res.y || z2 >= m_res.z) {
			/* Do an integer bounds test (may seem redundant - this is
			   to avoid a segfault, should a NaN/Inf ever find its way here..) */
			return 0;
		}

		const Float fx = p.x-x1, fy = p.y-y1, fz = p.z-z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f-fz;

		const Float
			d000 = m_data[(z1*m_res.y + y1)*m_res.x + x1],
			d001 = m_data[(z1*m_res.y + y1)*m_res.x + x2],
			d010 = m_data[(z1*m_res.y + y2)*m_res.x + x1],
			d011 = m_data[(z1*m_res.y + y2)*m_res.x + x2],
			d100 = m_data[(z2*m_res.y + y1)*m_res.x + x1],
			d101 = m_data[(z2*m_res.y + y1)*m_res.x + x2],
			d110 = m_data[(z2*m_res.y + y2)*m_res.x + x1],
			d111 = m_data[(z2*m_res.y + y2)*m_res.x + x2];

		return ((d000*_fx + d001*fx)*_fy +
				(d010*_fx + d011*fx)*fy)*_fz +
				((d100*_fx + d101*fx)*_fy +
				(d110*_fx + d111*fx)*fy)*fz;
	}


	Spectrum lookupSpectrum(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		if (p.x < 0 || p.y < 0 || p.z < 0)
			return Spectrum(0.0f);

		const int x1 = (int) p.x, y1 = (int) p.y, z1 = (int) p.z,
				x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x || 
			y2 >= m_res.y || z2 >= m_res.z) {
			/* Do an integer bounds test (may seem redundant - this is
			   to avoid a segfault, should a NaN/Inf ever find its way here..) */
			return Spectrum(0.0f);
		}

		const Float fx = p.x-x1, fy = p.y-y1, fz = p.z-z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f-fz;

		const float3 *spectrumData = (float3 *) m_data;

		const float3
			&d000 = spectrumData[(z1*m_res.y + y1)*m_res.x + x1],
			&d001 = spectrumData[(z1*m_res.y + y1)*m_res.x + x2],
			&d010 = spectrumData[(z1*m_res.y + y2)*m_res.x + x1],
			&d011 = spectrumData[(z1*m_res.y + y2)*m_res.x + x2],
			&d100 = spectrumData[(z2*m_res.y + y1)*m_res.x + x1],
			&d101 = spectrumData[(z2*m_res.y + y1)*m_res.x + x2],
			&d110 = spectrumData[(z2*m_res.y + y2)*m_res.x + x1],
			&d111 = spectrumData[(z2*m_res.y + y2)*m_res.x + x2];

		return (((d000*_fx + d001*fx)*_fy +
				(d010*_fx + d011*fx)*fy)*_fz +
				((d100*_fx + d101*fx)*_fy +
				(d110*_fx + d111*fx)*fy)*fz).toSpectrum();
	}

	Vector lookupVector(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		if (p.x < 0 || p.y < 0 || p.z < 0)
			return Vector(0.0f);

		const int x1 = (int) p.x, y1 = (int) p.y, z1 = (int) p.z,
				x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x || 
			y2 >= m_res.y || z2 >= m_res.z) {
			/* Do an integer bounds test (may seem redundant - this is
			   to avoid a segfault, should a NaN/Inf ever find its way here..) */
			return Vector(0.0f);
		}

		const Float fx = p.x-x1, fy = p.y-y1, fz = p.z-z1;
		const float3 *vectorData = (float3 *) m_data;
#if 1
		/* Nearest neighbor */
		return m_volumeToWorld(vectorData[
			(((fz < .5) ? z1 : z2)  * m_res.y +
			((fy < .5) ? y1 : y2)) * m_res.x +
			((fx < .5) ? x1 : x2)].toVector());
#else
		Float _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f-fz;

		const float3
			&d000 = vectorData[(z1*m_res.y + y1)*m_res.x + x1],
			&d001 = vectorData[(z1*m_res.y + y1)*m_res.x + x2],
			&d010 = vectorData[(z1*m_res.y + y2)*m_res.x + x1],
			&d011 = vectorData[(z1*m_res.y + y2)*m_res.x + x2],
			&d100 = vectorData[(z2*m_res.y + y1)*m_res.x + x1],
			&d101 = vectorData[(z2*m_res.y + y1)*m_res.x + x2],
			&d110 = vectorData[(z2*m_res.y + y2)*m_res.x + x1],
			&d111 = vectorData[(z2*m_res.y + y2)*m_res.x + x2];

		return m_volumeToWorld((((d000*_fx + d001*fx)*_fy +
				(d010*_fx + d011*fx)*fy)*_fz +
				((d100*_fx + d101*fx)*_fy +
				(d110*_fx + d111*fx)*fy)*fz).toVector());
#endif
	}
	
	bool supportsFloatLookups() const {
		return m_channels == 1;
	}
	
	bool supportsSpectrumLookups() const {
		return m_channels == 3;
	}
	
	bool supportsVectorLookups() const {
		return m_channels == 3;
	}

	Float getStepSize() const {
		return m_stepSize;
	}

	MTS_DECLARE_CLASS()
protected:
	fs::path m_filename;
	float *m_data;
	bool m_fromStream, m_sendData;
	Vector3i m_res;
	int m_channels;
	Transform m_worldToGrid;
	Transform m_worldToVolume;
	Transform m_volumeToWorld;
	Float m_stepSize;
	bool m_overrideAABB;
	AABB m_originalAABB;
#if defined(__LINUX__) || defined(__OSX__)
	size_t m_mmapSize;
	void *m_mmapPtr;
#elif defined(WIN32)
	HANDLE m_file;
	HANDLE m_fileMapping;
	void *m_mmapPtr;
#endif
};

MTS_IMPLEMENT_CLASS_S(GridDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(GridDataSource, "Grid data source");
MTS_NAMESPACE_END
