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

