#pragma once

#include "geometry_lib.h"

/**
 * @brief The Polygon struct
 * Convex polygon
 */
class Polygon
{
public:
	Point3f arr[MAX_VERTEX_NUM];
	size_t size = 0;

	Polygon();
	explicit Polygon(int size);
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
	size_t size = 0;

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
