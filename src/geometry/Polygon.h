#pragma once

#include "Point.h"
#include <iostream>

#define MIN_VERTEX_NUM 3		///< minimum number of vertices in polygon
#define MAX_VERTEX_NUM 64		///< maximum number of vertices in polygon

#define MAX_POLYGON_NUM 512		///< maximum number of polygons in array of polygons

/**
 * @brief The Polygon struct
 * Convex polygon
 */
class Polygon
{
public:
	Point3f arr[MAX_VERTEX_NUM];
	int nVertices = 0;

	Polygon();
	explicit Polygon(int nVertices);
	Polygon(const Polygon &other);
	Polygon(Polygon &&other);

	Polygon & operator = (const Polygon &other);
	Polygon & operator = (Polygon &&other);
	friend std::ostream & operator << (std::ostream &os, const Polygon &beam);

	double Area() const;
	Point3f Center() const;
	Point3f Normal() const;
};

class PolygonArray
{
public:
	Polygon arr[MAX_POLYGON_NUM];
	int size = 0;

	void Push(const Polygon &p)
	{
		arr[size++] = p;
	}

	Polygon &Pop()
	{
		return arr[--size];
	}

	void Clear()
	{
		size = 0;
	}
};
