/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/volume.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>
#include <mitsuba/core/mmap.h>


// Uncomment to enable nearest-neighbor direction interpolation
// #define VINTERP_NEAREST_NEIGHBOR

// Uncomment to enable linear direction interpolation (usually a bad idea)
//#define VINTERP_LINEAR

// Uncomment to enable structure tensor-based direction interpolation (best, but slowest)
#define VINTERP_STRUCTURE_TENSOR

// Number of power iteration steps used to find the dominant direction
#define POWER_ITERATION_STEPS 5

MTS_NAMESPACE_BEGIN

/**
 * \brief This class implements access to volume data stored on a 
 * 3D grid using a simple binary exchange format.
 *
 * The format uses a little endian encoding and is specified as follows:
 *
 * Bytes 1-3   :  ASCII Bytes 'V', 'O', and 'L' 
 * Byte  4     :  Version identifier (currently 3)
 * Bytes 5-8   :  Encoding identifier using a double word
 *
 *                    1 => Dense float32-based representation
 *
 *                    2 => Dense float16-based representation 
 *                         NOT YET SUPPORTED BY THIS IMPLEMENTATION
 *
 *                    3 => Dense uint8-based representation
 *                         The range 0..255 will be mapped to 0..1.
 *
 *                    4 => Dense quantized directions
 *                         The directions are stored in spherical coordinates,
 *                         with a total storage cost of 16 bit per entry.
 *
 * Bytes 9-12  :  Number of cells along the X axis (double word)
 * Bytes 13-16 :  Number of cells along the Y axis (double word)
 * Bytes 17-20 :  Number of cells along the Z axis (double word)
 * Bytes 21-24 :  Number of channels (double word, supported values: 1 or 3)
 * Bytes 25-48 :  Axis-aligned bounding box of the data stored in single
 *                precision (order: xmin, ymin, zmin, xmax, ymax, zmax)
 * Bytes 49-*  :  Binary data of the volume stored in the specified encoding.
 *                The data are ordered so that the following C-style indexing
 *                operation makes sense after the file has been mapped into memory:
 *                   "data[((zpos*yres + ypos)*xres + xpos)*channels + chan]"
 *                where (xpos, ypos, zpos, chan) denotes the lookup location.
 *
 * Note that Mitsuba expects that entries in direction volumes are either
 * zero or valid unit vectors.
 */
class GridDataSource : public VolumeDataSource {
public:
	enum EVolumeType {
		EFloat32 = 1,
		EFloat16 = 2,
		EUInt8 = 3,
		EQuantizedDirections = 4
	};

	GridDataSource(const Properties &props) 
		: VolumeDataSource(props) {
		m_volumeToWorld = props.getTransform("toWorld", Transform());

		if (props.hasProperty("min") && props.hasProperty("max")) {
			/* Optionally allow to use an AABB other than 
			   the one specified by the grid file */
			m_dataAABB.min = props.getPoint("min");
			m_dataAABB.max = props.getPoint("max");
		}

		/**
		 * When 'sendData' is set to false, only the filename 
		 * is transmitted. A following unserialization of the 
		 * stream causes the implementation to then look for 
		 * the file (which had better exist if unserialization 
		 * occurs on a remote machine). 
		 */
		 m_sendData = props.getBoolean("sendData", false);

		 loadFromFile(props.getString("filename"));
	}

	GridDataSource(Stream *stream, InstanceManager *manager) 
			: VolumeDataSource(stream, manager) {
		m_volumeToWorld = Transform(stream);
		m_dataAABB = AABB(stream);
		m_sendData = stream->readBool();
		if (m_sendData) { 
			m_volumeType = (EVolumeType) stream->readInt();
			m_res = Vector3i(stream);
			m_channels = stream->readInt();
			m_filename = stream->readString();
			size_t volumeSize = getVolumeSize();
			m_data = new uint8_t[volumeSize];
			stream->read(m_data, volumeSize);
		} else {
			fs::path filename = stream->readString();
			loadFromFile(filename);
		}
		configure();
	}

	virtual ~GridDataSource() {
		if (!m_mmap)
			delete[] m_data;
	}

	size_t getVolumeSize() const {
		size_t nEntries = (size_t) m_res.x 
			* (size_t) m_res.y * (size_t) m_res.z;
		switch (m_volumeType) {
			case EFloat32: return 4 * nEntries * m_channels;
			case EFloat16: return 2 * nEntries * m_channels;
			case EUInt8:   return 1 * nEntries * m_channels;
			case EQuantizedDirections:  return 2 * nEntries;
			default:
				Log(EError, "Unknown volume format!");
				return 0;
		}
	}

	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);

		m_volumeToWorld.serialize(stream);
		m_dataAABB.serialize(stream);
		stream->writeBool(m_sendData);

		if (m_sendData) {
			stream->writeInt(m_volumeType);
			m_res.serialize(stream);
			stream->writeInt(m_channels);
			stream->writeString(m_filename.file_string());
			stream->write(m_data, getVolumeSize());
		} else {
			stream->writeString(m_filename.file_string());
		}
	}

	void configure() {
		Vector extents(m_dataAABB.getExtents());
		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(
				(m_res[0] - 1) / extents[0],
				(m_res[1] - 1) / extents[1],
				(m_res[2] - 1) / extents[2])
			) * Transform::translate(-Vector(m_dataAABB.min)) * m_worldToVolume;
		m_stepSize = std::numeric_limits<Float>::infinity();
		for (int i=0; i<3; ++i)
			m_stepSize = 0.5f * std::min(m_stepSize, extents[i] / (Float) (m_res[i]-1));
		m_aabb.reset();
		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_volumeToWorld(m_dataAABB.getCorner(i)));

		/* Precompute cosine and sine lookup tables */
		for (int i=0; i<255; i++) {
			Float angle = (float) i * ((float) M_PI / 255.0f);
			m_cosPhi[i] = std::cos(2.0f * angle);
			m_sinPhi[i] = std::sin(2.0f * angle);
			m_cosTheta[i] = std::cos(angle);
			m_sinTheta[i] = std::sin(angle);
			m_densityMap[i] = i/255.0f;
		}
		m_cosPhi[255] = m_sinPhi[255] = 0;
		m_cosTheta[255] = m_sinTheta[255] = 0;
		m_densityMap[255] = 1.0f;
	}

	void loadFromFile(const fs::path &filename) {
		m_filename = filename;
		fs::path resolved = Thread::getThread()->getFileResolver()->resolve(filename);
		m_mmap = new MemoryMappedFile(resolved);
		ref<MemoryStream> stream = new MemoryStream(m_mmap->getData(), m_mmap->getSize());
		stream->setByteOrder(Stream::ELittleEndian);

		char header[3];
		stream->read(header, 3);
		if (header[0] != 'V' || header[1] != 'O' || header[2] != 'L')
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect header identifier)");
		uint8_t version;
		stream->read(&version, 1);
		if (version != 3)
			Log(EError, "Encountered an invalid volume data file "
				"(incorrect file version)");
		int type = stream->readInt();

		int xres = stream->readInt(),
			yres = stream->readInt(),
			zres = stream->readInt();
		m_res = Vector3i(xres, yres, zres);
		m_channels = stream->readInt();

		switch (type) {
			case EFloat32:
				if (m_channels != 1 && m_channels != 3)
					Log(EError, "Encountered an unsupported float32 volume data "
						"file (%i channels, only 1 and 3 are supported)",
						m_channels);
				break;
			case EFloat16:
				Log(EError, "Error: float16 volumes are not yet supported!");
			case EUInt8:
				if (m_channels != 1 && m_channels != 3)
					Log(EError, "Encountered an unsupported uint8 volume data "
						"file (%i channels, only 1 and 3 are supported)", m_channels);
				break;
			case EQuantizedDirections:
				if (m_channels != 3)
					Log(EError, "Encountered an unsupported quantized direction "
							"volume data file (%i channels, only 3 are supported)",
							m_channels);
				break;
			default:
				Log(EError, "Encountered a volume data file of unknown type!");
		}

		m_volumeType = (EVolumeType) type;

		if (!m_dataAABB.isValid()) {
			Float xmin = stream->readSingle(),
				  ymin = stream->readSingle(),
				  zmin = stream->readSingle();
			Float xmax = stream->readSingle(),
				  ymax = stream->readSingle(),
				  zmax = stream->readSingle();
			m_dataAABB = AABB(Point(xmin, ymin, zmin), Point(xmax, ymax, zmax));
		}

		Log(EDebug, "Mapped \"%s\" into memory: %ix%ix%i (%i channels), %s, %s", 
			resolved.filename().c_str(), m_res.x, m_res.y, m_res.z, m_channels,
			memString(m_mmap->getSize()).c_str(), m_dataAABB.toString().c_str());
		m_data = (uint8_t *) (((float *) m_mmap->getData()) + 12);
	}

	/**
	 * This is needed since Mitsuba might be 
	 * compiled with either single/double precision
	 */
	struct float3 {
		float value[3];

		inline float3() { }

		inline float3(float a, float b, float c) {
			value[0] = a; value[1] = b; value[2] = c;
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
	
		float operator[](int i) const {
			return value[i];
		}

		inline Matrix3x3 tensor() const {
			return Matrix3x3(
				value[0]*value[0], value[0]*value[1], value[0]*value[2],
				value[1]*value[0], value[1]*value[1], value[1]*value[2],
				value[2]*value[0], value[2]*value[1], value[2]*value[2]
			);
		}
	};

	Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = floorToInt(p.x),
			  y1 = floorToInt(p.y),
			  z1 = floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z) 
			return 0;

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
				const float *floatData = (float *) m_data;
				const Float
					d000 = floatData[(z1*m_res.y + y1)*m_res.x + x1],
					d001 = floatData[(z1*m_res.y + y1)*m_res.x + x2],
					d010 = floatData[(z1*m_res.y + y2)*m_res.x + x1],
					d011 = floatData[(z1*m_res.y + y2)*m_res.x + x2],
					d100 = floatData[(z2*m_res.y + y1)*m_res.x + x1],
					d101 = floatData[(z2*m_res.y + y1)*m_res.x + x2],
					d110 = floatData[(z2*m_res.y + y2)*m_res.x + x1],
					d111 = floatData[(z2*m_res.y + y2)*m_res.x + x2];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			case EUInt8: {
				const Float
					d000 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x1]],
					d001 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x2]],
					d010 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x1]],
					d011 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x2]],
					d100 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x1]],
					d101 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x2]],
					d110 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x1]],
					d111 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x2]];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			default:
				return 0.0f;
		}
	}

	Spectrum lookupSpectrum(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = floorToInt(p.x),
			  y1 = floorToInt(p.y),
			  z1 = floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z) 
			return Spectrum(0.0f);

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
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
			case EUInt8: {
				const float3
					d000 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+2]]),
					d001 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+2]]),
					d010 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+2]]),
					d011 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+2]]),
					d100 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+2]]),
					d101 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+2]]),
					d110 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+2]]),
					d111 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+2]]);

				return (((d000*_fx + d001*fx)*_fy +
						 (d010*_fx + d011*fx)*fy)*_fz +
						((d100*_fx + d101*fx)*_fy +
						 (d110*_fx + d111*fx)*fy)*fz).toSpectrum();

				}
			default: return Spectrum(0.0f);
		}
	}

	Vector lookupVector(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = floorToInt(p.x),
			  y1 = floorToInt(p.y),
			  z1 = floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z) 
			return Vector(0.0f);

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1;
		Vector value;

		#if defined(VINTERP_NEAREST_NEIGHBOR)
			/* Nearest neighbor */
			switch (m_volumeType) {
				case EFloat32: {
					const float3 *vectorData = (float3 *) m_data;
					value = vectorData[
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2)].toVector();
					}
					break;
				case EQuantizedDirections: {
					value = lookupQuantizedDirection(
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2)).toVector();
					}
					break;
				default:
					return Vector(0.0f);
			}
		#else
			Float _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;
			float3 d[8];

			switch (m_volumeType) {
				case EFloat32: {
						const float3 *vectorData = (float3 *) m_data;
						d[0] = vectorData[(z1*m_res.y + y1)*m_res.x + x1];
						d[1] = vectorData[(z1*m_res.y + y1)*m_res.x + x2];
						d[2] = vectorData[(z1*m_res.y + y2)*m_res.x + x1];
						d[3] = vectorData[(z1*m_res.y + y2)*m_res.x + x2];
						d[4] = vectorData[(z2*m_res.y + y1)*m_res.x + x1];
						d[5] = vectorData[(z2*m_res.y + y1)*m_res.x + x2];
						d[6] = vectorData[(z2*m_res.y + y2)*m_res.x + x1];
						d[7] = vectorData[(z2*m_res.y + y2)*m_res.x + x2];
					}
					break;
				case EQuantizedDirections: {
						d[0] = lookupQuantizedDirection((z1*m_res.y + y1)*m_res.x + x1);
						d[1] = lookupQuantizedDirection((z1*m_res.y + y1)*m_res.x + x2);
						d[2] = lookupQuantizedDirection((z1*m_res.y + y2)*m_res.x + x1);
						d[3] = lookupQuantizedDirection((z1*m_res.y + y2)*m_res.x + x2);
						d[4] = lookupQuantizedDirection((z2*m_res.y + y1)*m_res.x + x1);
						d[5] = lookupQuantizedDirection((z2*m_res.y + y1)*m_res.x + x2);
						d[6] = lookupQuantizedDirection((z2*m_res.y + y2)*m_res.x + x1);
						d[7] = lookupQuantizedDirection((z2*m_res.y + y2)*m_res.x + x2);
					}
					break;
				default:
					return Vector(0.0f);
			}

			#if defined(VINTERP_LINEAR)
				value = (((d[0]*_fx + d[1]*fx)*_fy +
						  (d[2]*_fx + d[3]*fx)* fy)*_fz +
						 ((d[4]*_fx + d[5]*fx)*_fy +
						  (d[6]*_fx + d[7]*fx)* fy)* fz).toVector();
			#elif defined(VINTERP_STRUCTURE_TENSOR)
				Matrix3x3 tensor(0.0f);
				for (int k=0; k<8; ++k) {
					Float factor = ((k & 1) ? fx : _fx) * ((k & 2) ? fy : _fy) 
						* ((k & 4) ? fz : _fz);
					for (int i=0; i<3; ++i) 
						for (int j=0; j<3; ++j) 
							tensor(i, j) += factor * d[k][i] * d[k][j];
				}
				if (tensor.isZero())
					return Vector(0.0f);
	
#if 0
				Float lambda[3];
				eig3_noniter(tensor, lambda);
				value = tensor.col(0);
#else

				/* Square the structure tensor for faster convergence */
				tensor *= tensor;

				const Float invSqrt3 = 0.577350269189626f;
				value = Vector(invSqrt3, invSqrt3, invSqrt3);

				/* Determine the dominant eigenvector using
				a few power iterations */
				for (int i=0; i<POWER_ITERATION_STEPS-1; ++i)
					value = normalize(tensor * value);
				value = tensor * value;
#endif

//				Float specularity = 1-lambda[1]/lambda[0];
			#else
				#error Need to choose a vector interpolation method!
			#endif
		#endif

		if (!value.isZero())
			return normalize(m_volumeToWorld(value));
		else
			return Vector(0.0f);
	}

	bool supportsFloatLookups() const { return m_channels == 1; }
	bool supportsSpectrumLookups() const { return m_channels == 3; }
	bool supportsVectorLookups() const { return m_channels == 3; }
	Float getStepSize() const { return m_stepSize; }

	MTS_DECLARE_CLASS()
protected: 
	inline float3 lookupQuantizedDirection(size_t index) const {
		uint8_t theta = m_data[2*index], phi = m_data[2*index+1];
		return float3(
			m_cosPhi[phi] * m_sinTheta[theta],
			m_sinPhi[phi] * m_sinTheta[theta],
			m_cosTheta[theta]
		);
	}

protected:
	fs::path m_filename;
	uint8_t *m_data;
	bool m_sendData;
	EVolumeType m_volumeType;
	Vector3i m_res;
	int m_channels;
	Transform m_worldToGrid;
	Transform m_worldToVolume;
	Transform m_volumeToWorld;
	Float m_stepSize;
	AABB m_dataAABB;
	ref<MemoryMappedFile> m_mmap;
	float m_cosTheta[256], m_sinTheta[256];
	float m_cosPhi[256], m_sinPhi[256];
	Float m_densityMap[256];
};

MTS_IMPLEMENT_CLASS_S(GridDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(GridDataSource, "Grid data source");
MTS_NAMESPACE_END
