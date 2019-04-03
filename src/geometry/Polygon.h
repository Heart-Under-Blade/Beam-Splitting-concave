#pragma once

#include "Point.h"
#include <iostream>

#define MIN_VERTEX_NUM 3		///< minimum number of vertices in polygon

#ifdef _DEBUG // DEB
#define MAX_VERTEX_NUM 256		///< maximum number of vertices in polygon
#else
/// REF OPT: уменьшить в 2 раза
#define MAX_VERTEX_NUM 64		///< maximum number of vertices in polygon
#endif

#define MAX_POLYGON_NUM 512		///< maximum number of polygons in array of polygons

/**
 * @brief Polygon consisted of 3 or more 3-coordinated vertices.
 * It will work correctly if the points are situated on the same plane.
 */
class Polygon
{
public:
	Point3f vertices[MAX_VERTEX_NUM];
	int nVertices = 0;

	Polygon();
	explicit Polygon(int nVertices);
	Polygon(const Polygon &other);
	Polygon(Polygon &&other);

	Polygon & operator = (const Polygon &other);
	Polygon & operator = (Polygon &&other);

	void AddVertex(const Point3f &v);
	void InsertVertex(int index, const Point3f &v);
	void RemoveVertex(int index);
	void Concat(const Polygon &other);
	void Clear();

	double Area() const;
	Point3f Center() const;
	Point3f Normal() const;

	friend std::ostream & operator << (std::ostream &os, const Polygon &beam);
};

class PolygonStack
{
public:
	Polygon polygons[MAX_POLYGON_NUM];
	int nPolygons = 0;

	void Push(const Polygon &p)
	{
		polygons[nPolygons++] = p;
	}

	Polygon &Pop()
	{
		return polygons[--nPolygons];
	}

	void Clear()
	{
		nPolygons = 0;
	}
};
