#ifndef HAIR_H_
#define HAIR_H_

#include <mitsuba/render/shape.h>

MTS_NAMESPACE_BEGIN


class Hair : public Shape {

	// The radius of the hair fibers (all fibers have constant radius)
	Float m_radius;

	// The vertices of the hair: [0...nVertices)
	std::vector<Point> m_vertices;

	// An indication of which vertices start a new fiber: [0...nVertices+1)
	std::vector<bool> m_startFiber;

	// A mapping of segment indices to vertex indices (needed only for construction): [0...nSegments)
	std::vector<int> m_segIndex;

public:

	Hair(const Properties &props);

	Hair(Stream *stream, InstanceManager *manager);

	void serialize(Stream *stream, InstanceManager *manager) const;

	bool isCompound() const {
		return true;
	}

	Shape *getElement(int index);

	Float radius() const { return m_radius; }

	Point vertex(int iv) const { return m_vertices[iv]; }
	void vertex(int iv, Point &p) const { p = m_vertices[iv]; }

	bool vertexStartsFiber(int iv) const { return m_startFiber[iv]; }


protected:

	void buildSegIndex();

	MTS_DECLARE_CLASS()
};


MTS_NAMESPACE_END


#endif /* HAIR_H_ */
