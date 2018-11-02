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
 * @brief Polygon consisted of 3-coordinate vertices
 */
class Polygon1
{
public:
	Point3f arr[MAX_VERTEX_NUM];
	int nVertices = 0;

	Polygon1();
	explicit Polygon1(int nVertices);
	Polygon1(const Polygon1 &other);
	Polygon1(Polygon1 &&other);

	void AddVertex(const Point3f &v);

	void InsertVertex(int index, const Point3f &v);
	void DeleteVertex(int index);

	Polygon1 & operator = (const Polygon1 &other);
	Polygon1 & operator = (Polygon1 &&other);
	friend std::ostream & operator << (std::ostream &os, const Polygon1 &beam);

	double Area() const;
	Point3f Center() const;
	Point3f Normal() const;
};

class PolygonArray
{
public:
	Polygon1 arr[MAX_POLYGON_NUM];
	int size = 0;

	void Push(const Polygon1 &p)
	{
		arr[size++] = p;
	}

	Polygon1 &Pop()
	{
		return arr[--size];
	}

	void Clear()
	{
		size = 0;
	}
};
