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

#if !defined(__RFILTER_H)
#define __RFILTER_H

#include <mitsuba/core/cobject.h>
#include <mitsuba/core/properties.h>

MTS_NAMESPACE_BEGIN

/**
 * \brief Abstract image reconstruction filter
 * \ingroup librender
 */
class MTS_EXPORT_RENDER ReconstructionFilter : public ConfigurableObject {
public:
	/// Return the filter's width
	inline const Vector2 &getFilterSize() const {
		return m_size;
	}

	/// Evaluate the filter function
	virtual Float evaluate(Float x, Float y) const = 0;

	/// Sample from the filter. Returns a weight
	virtual Float sample(const Point2 &sample, Float &x, Float &y) const;

	/// Return the properties of this reconstruction filter
	inline const Properties &getProperties() const { return m_properties; }

	MTS_DECLARE_CLASS()
protected:
    /// Create a new reconstruction filter
    ReconstructionFilter(const Properties &props);

	/// Unserialize a filter
	ReconstructionFilter(Stream *stream, InstanceManager *manager);

	/// Virtual destructor
	virtual ~ReconstructionFilter();
protected:
	Vector2 m_size;
	Properties m_properties;
};

/// Filter resolution for \ref TabulatedFilter
#define FILTER_RESOLUTION 15

/**
 * \brief Tabulates expensive-to-evaluate filters so that they become
 * simple array lookups. 
 *
 * Only works for symmetric filters!
 * \ingroup librender
 */
class TabulatedFilter : public Object {
public:
	/// Tabulate a reconstruction filter
    TabulatedFilter(const ReconstructionFilter *filter);

	/// Unserialize a filter
	TabulatedFilter(Stream *stream);

	/// Return the name of the tabulated filter
	inline const std::string &getName() const {
		return m_name;
	}

	/// Return the filter's width
	inline const Vector2 &getFilterSize() const {
		return m_size;
	}
	
	/// Return the size factor of the underlying discretization
	inline const Vector2 &getSizeFactor() const {
		return m_factor;
	}

	/// Evaluate the filter function
	inline Float evaluate(Float x, Float y) const {
		int xPos = (int) (m_factor.x * std::abs(x));
		int yPos = (int) (m_factor.y * std::abs(y));
		if (xPos >= FILTER_RESOLUTION || yPos >= FILTER_RESOLUTION)
			return 0.0f;
		else
			return m_values[yPos][xPos];
	}

	/// Perform a lookup into the underlying array discretization
	inline Float lookup(int x, int y) const {
		return m_values[y][x];
	}

	/// Serialize a filter
	void serialize(Stream *stream) const;
    
	/// Return a string representation
	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~TabulatedFilter();
protected:
	Vector2 m_size, m_factor;
	std::string m_name;
	Float m_values[FILTER_RESOLUTION+1][FILTER_RESOLUTION+1];
};

MTS_NAMESPACE_END

#endif /* __RFILTER_H */
