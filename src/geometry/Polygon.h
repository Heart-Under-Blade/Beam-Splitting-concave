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
	size_t nVertices = 0;

	Polygon();
	explicit Polygon(int nVertices);
	Polygon(const Polygon &other);
	Polygon(Polygon &&other);

	void AddVertex(const Point3f &v);

	void Concat(const Polygon &other)
	{
		for (int i = 0; i < other.nVertices; ++i)
		{
			AddVertex(other.arr[i]);
		}
	}

	void InsertVertex(int index, const Point3f &v)
	{
		++nVertices;

		for (int i = nVertices-1; i > index; --i)
		{
			arr[i] = arr[i-1];
		}

		arr[index] = v;
	}

	void DeleteVertex(int index)
	{
		for (int i = nVertices-1; i > index; --i)
		{
			arr[i-1] = arr[i];
		}

		--nVertices;
	}

	Polygon & operator = (const Polygon &other);
	Polygon & operator = (Polygon &&other);
	friend std::ostream & operator << (std::ostream &os, const Polygon &beam);

	void Clear()
	{
		nVertices = 0;
	}

	double Area() const;
	Point3f Center() const;
	Point3f Normal() const;
};

class BigPolygon
{
public:
	Point3f arr[16384];
	int nVertices = 0;

	BigPolygon() {}

	void AddVertex(const Point3f &v)
	{
		arr[nVertices++] = v;
	}

	void Concat(const Polygon &other)
	{
		for (int i = 0; i < other.nVertices; ++i)
		{
			AddVertex(other.arr[i]);
		}
	}
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
