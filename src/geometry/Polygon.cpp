#include "Polygon.h"

Polygon::Polygon()
{

}

Polygon::Polygon(int size) : size(size) {}

Polygon::Polygon(const Polygon &other)
{
	size = other.size;

	for (int i = 0; i < other.size; ++i)
	{
		arr[i] = other.arr[i];
	}
}

Polygon::Polygon(Polygon &&other)
{
	size = other.size;

	for (int i = 0; i < size; ++i)
	{
		arr[i] = other.arr[i];
	}

	other.size = 0;
}

Polygon &Polygon::operator =(const Polygon &other)
{
	if (this != &other)
	{
		size = other.size;

		for (int i = 0; i < size; ++i)
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
		size = other.size;

		for (int i = 0; i < size; ++i)
		{
			arr[i] = other.arr[i];
		}

		other.size = 0;
	}

	return *this;
}

std::ostream &operator <<(std::ostream &os, const Polygon &beam)
{
	using namespace std;

	os << "polygon: {" << endl;

	for (int i = 0; i < beam.size; ++i)
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

	for (int i = 2; i < size; ++i)
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

	for (int i = 0; i < size; ++i)
	{
		p = p + arr[i];
	}

	return p/size;
}

Point3f Polygon::Normal() const
{
	Point3f normal;

	Point3f p1 = arr[1] - arr[0];
	Point3f p2 = arr[2] - arr[0];
	CrossProduct(p1, p2, normal);

	Normalize(normal);
	return normal;
}
