#include "Polygon.h"
#include <math.h>

#define EPS_NORMAL 0.1

Polygon1::Polygon1()
{

}

Polygon1::Polygon1(int size) : nVertices(size) {}

Polygon1::Polygon1(const Polygon1 &other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < other.nVertices; ++i)
	{
		arr[i] = other.arr[i];
	}
}

Polygon1::Polygon1(Polygon1 &&other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < nVertices; ++i)
	{
		arr[i] = other.arr[i];
	}

	other.nVertices = 0;
}

void Polygon1::AddVertex(const Point3f &v)
{
	arr[nVertices++] = v;
}

void Polygon1::InsertVertex(int index, const Point3f &v)
{
	++nVertices;

	for (int i = nVertices-1; i > index; --i)
	{
		arr[i] = arr[i-1];
	}

	arr[index] = v;
}

void Polygon1::DeleteVertex(int index)
{
	for (int i = nVertices-1; i > index; --i)
	{
		arr[i-1] = arr[i];
	}

	--nVertices;
}

Polygon1 &Polygon1::operator =(const Polygon1 &other)
{
	if (this != &other)
	{
		nVertices = other.nVertices;

		for (int i = 0; i < nVertices; ++i)
		{
			arr[i] = other.arr[i];
		}
	}

	return *this;
}

Polygon1 &Polygon1::operator = (Polygon1 &&other)
{
	if (this != &other)
	{
		nVertices = other.nVertices;

		for (int i = 0; i < nVertices; ++i)
		{
			arr[i] = other.arr[i];
		}

		other.nVertices = 0;
	}

	return *this;
}

std::ostream &operator <<(std::ostream &os, const Polygon1 &beam)
{
	for (int i = 0; i < beam.nVertices; ++i)
	{
		os /*<< "\t"*/
		   << beam.arr[i].coordinates[0] << " "
		   << beam.arr[i].coordinates[1] << " "
		   << beam.arr[i].coordinates[2] << " " << std::endl;
	}

	os << beam.arr[0].coordinates[0] << " "
	   << beam.arr[0].coordinates[1] << " "
	   << beam.arr[0].coordinates[2] << " " << std::endl;

	return os;
}

double Polygon1::Area() const
{
	double square = 0;
	const Point3f &basePoint = arr[0];
	Point3f p1 = arr[1] - basePoint;

	for (int i = 2; i < nVertices; ++i)
	{
		Point3f p2 = arr[i] - basePoint;
		Point3f res;
		Point3f::CrossProduct(p1, p2, res);
		square += Point3f::Length(res);
		p1 = p2;
	}

	return square/2.0;
}

Point3f Polygon1::Center() const
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
		p = p + arr[i];
	}

	return p/nVertices;
}

Point3f Polygon1::Normal() const
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
		Point3f p1 = arr[first] - arr[start];
		Point3f p2 = arr[next] - arr[start];
		Point3f::CrossProduct(p1, p2, normal);
		Point3f::Normalize(normal);
		++count;
	}
	while (isnan(normal.coordinates[0]) && count < nVertices);

	return normal;
}
