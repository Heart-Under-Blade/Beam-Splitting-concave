#include "Polygon.h"
#include <math.h>

#define EPS_NORMAL 0.1

Polygon::Polygon()
{
}

Polygon::Polygon(int size) : nVertices(size) {}

Polygon::Polygon(const Polygon &other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < other.nVertices; ++i)
	{
		vertices[i] = other.vertices[i];
	}
}

Polygon::Polygon(Polygon &&other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < nVertices; ++i)
	{
		vertices[i] = other.vertices[i];
	}

	other.nVertices = 0;
}

void Polygon::AddVertex(const Point3f &v)
{
	vertices[nVertices++] = v;
}

void Polygon::Concat(const Polygon &other)
{
	for (int i = 0; i < other.nVertices; ++i)
	{
		AddVertex(other.vertices[i]);
	}
}

void Polygon::InsertVertex(int index, const Point3f &v)
{
	++nVertices;

	for (int i = nVertices-1; i > index; --i)
	{
		vertices[i] = vertices[i-1];
	}

	vertices[index] = v;
}

void Polygon::RemoveVertex(int index)
{
	for (int i = index+1; i < nVertices; ++i)
	{
		vertices[i-1] = vertices[i];
	}

	--nVertices;
}

Polygon &Polygon::operator =(const Polygon &other)
{
	if (this != &other)
	{
		nVertices = other.nVertices;

		for (int i = 0; i < nVertices; ++i)
		{
			vertices[i] = other.vertices[i];
		}
	}

	return *this;
}

Polygon &Polygon::operator = (Polygon &&other)
{
	if (this != &other)
	{
		nVertices = other.nVertices;

		for (int i = 0; i < nVertices; ++i)
		{
			vertices[i] = other.vertices[i];
		}

		other.nVertices = 0;
	}

	return *this;
}

std::ostream &operator <<(std::ostream &os, const Polygon &p)
{
	for (int i = 0; i < p.nVertices; ++i)
	{
		os << p.vertices[i] << std::endl;
	}

	return os;
}

double Polygon::Area() const
{
	double square = 0;
	const Point3f &basePoint = vertices[0];
	Point3f p1 = vertices[1] - basePoint;

	for (int i = 2; i < nVertices; ++i)
	{
		Point3f p2 = vertices[i] - basePoint;
		Point3f res;
		Point3f::CrossProduct(p1, p2, res);
		square += Point3f::Length(res);
		p1 = p2;
	}

	return square/2.0;
}

Point3f Polygon::Center() const
{
#ifdef _DEBUG // DEB
	if (nVertices == 0)
	{
		std::cerr << "ERROR! Polygon is empty." << std::endl
				  << __FILE__ << ", " << __FUNCTION__ << std::endl;
		throw std::exception();
	}
#endif
	Point3f p(0, 0, 0);

	for (int i = 0; i < nVertices; ++i)
	{
		p = p + vertices[i];
	}

	return p/nVertices;
}

Point3f Polygon::Normal() const
{
	Point3f normal;

	int count = 0;
	int start, first, next;
	start = nVertices-1;

	do
	{
		start = (start + 1 != nVertices) ? start + 1 : 0;
		first = (start + 1 != nVertices) ? start + 1 : 0;
		next  = (first + 1 != nVertices) ? first + 1 : 0;
		Point3f p1 = vertices[first] - vertices[start];
		Point3f p2 = vertices[next] - vertices[start];
		Point3f::CrossProduct(p1, p2, normal);
		Point3f::Normalize(normal);
		++count;
	}
	while (isnan(normal.coordinates[0]) && count < nVertices);

	return normal;
}

void Polygon::Clear()
{
	nVertices = 0;
}
