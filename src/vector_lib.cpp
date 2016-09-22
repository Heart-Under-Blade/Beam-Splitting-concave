#include <math.h>
#include "types.h"

double DotProduct(const Point3f &v1, const Point3f &v2)
{
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z;
}

double Norm(const Point3f &point)
{
	return point.X*point.X + point.Y*point.Y + point.Z*point.Z;
}

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res)
{
	res.X = v1.Y*v2.Z - v1.Z*v2.Y;
	res.Y = v1.Z*v2.X - v1.X*v2.Z;
	res.Z = v1.X*v2.Y - v1.Y*v2.X;
}

void Normalize(Point3f &v)
{
	double lenght = sqrt(Norm(v));
	v.X /= lenght;
	v.Y /= lenght;
	v.Z /= lenght;
}

Point3f NormalToFacet(const Point3f *facet)
{
	Point3f normal;

	Point3f p1 = facet[1] - facet[0];
	Point3f p2 = facet[2] - facet[0];
	CrossProduct(p2, p1, normal);

	Normalize(normal);
	return normal;
}

void CopyPoints(Point3f *points, Point3f *result, int size)
{
	for (int i = 0; i <= size; ++i)
	{
		result[i].X = points[i].X;
		result[i].Y = points[i].Y;
		result[i].Z = points[i].Z;
	}
}
