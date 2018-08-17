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

struct PolygonArray
{
	Polygon arr[MAX_POLYGON_NUM];
	int size = 0;
};
