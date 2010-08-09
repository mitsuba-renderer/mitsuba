#if !defined(__IMAGEPROC_WU_H)
#define __IMAGEPROC_WU_H

#include <mitsuba/core/sched.h>

MTS_NAMESPACE_BEGIN

/**
 * Rectangular image work unit. Used by the <tt>BlockedImageProcess</tt>.
 */
class MTS_EXPORT_RENDER RectangularWorkUnit : public WorkUnit {
public:
	inline RectangularWorkUnit() { }

	/* WorkUnit implementation */
	void set(const WorkUnit *wu);
	void load(Stream *stream);
	void save(Stream *stream) const;

	inline const Point2i &getOffset() const { return m_offset; }
	inline const Vector2i &getSize() const { return m_size; }

	inline void setOffset(const Point2i &offset) { m_offset = offset; }
	inline void setSize(const Vector2i &size) { m_size = size; }

	std::string toString() const;

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~RectangularWorkUnit() { }
private:
	Point2i m_offset;
	Vector2i m_size;
};

MTS_NAMESPACE_END

#endif /* __IMAGEPROC_WU_H */
