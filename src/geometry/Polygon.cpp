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
		arr[i] = other.arr[i];
	}
}

Polygon::Polygon(Polygon &&other)
{
	nVertices = other.nVertices;

	for (int i = 0; i < nVertices; ++i)
	{
		arr[i] = other.arr[i];
	}

	other.nVertices = 0;
}

void Polygon::AddVertex(const Point3f &v)
{
	arr[nVertices++] = v;
}

Polygon &Polygon::operator =(const Polygon &other)
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

Polygon &Polygon::operator = (Polygon &&other)
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

std::ostream &operator <<(std::ostream &os, const Polygon &beam)
{
	using namespace std;

	os << "polygon: {" << endl;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		os << "\t" << i << ": "
		   << beam.arr[i].cx << ", "
		   << beam.arr[i].cy << ", "
		   << beam.arr[i].cz << ", "
		   << beam.arr[i].d_param << endl;
	}

	os << "}" << endl << endl;
	return os;
}

double Polygon::Area() const
{
	double square = 0;
	const Point3f &basePoint = arr[0];
	Point3f p1 = arr[1] - basePoint;

	for (int i = 2; i < nVertices; ++i)
	{
		Point3f p2 = arr[i] - basePoint;
		Point3f res;
		CrossProduct(p1, p2, res);
		square += Length(res);
		p1 = p2;
	}

	return square/2.0;
}

Point3f Polygon::Center() const
{
	Point3f p(0, 0, 0);

	for (int i = 0; i < nVertices; ++i)
	{
		p = p + arr[i];
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
		Point3f p1 = arr[first] - arr[start];
		Point3f p2 = arr[next] - arr[start];
		CrossProduct(p1, p2, normal);
		Normalize(normal);
		++count;
	}
	while (isnan(normal.cx) && count < nVertices);

	return normal;
}
