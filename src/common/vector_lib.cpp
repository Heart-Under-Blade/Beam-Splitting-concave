#include <math.h>
#include "types.h"

double DotProduct(const Point3f &v1, const Point3f &v2)
{
	return v1.cx*v2.cx + v1.cy*v2.cy + v1.cz*v2.cz;
}

double DotProductD(const Point3d &v1, const Point3d &v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

double Norm(const Point3f &point)
{
	return point.cx*point.cx + point.cy*point.cy + point.cz*point.cz;
}

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res)
{
	res.cx = v1.cy*v2.cz - v1.cz*v2.cy;
	res.cy = v1.cz*v2.cx - v1.cx*v2.cz;
	res.cz = v1.cx*v2.cy - v1.cy*v2.cx;
}

double Length(const Point3f &v)
{
	return sqrt(Norm(v));
}

void Normalize(Point3f &v) // OPT: всю функцию попробовать на интринсиках
{
	double lenght = sqrt(Norm(v)); // OPT: расписать
	v.cx /= lenght;
	v.cy /= lenght;
	v.cz /= lenght;
}

Point3f NormalToFacet(const Point3f *facet)
{
	Point3f normal;

	Point3f p1 = facet[1] - facet[0];
	Point3f p2 = facet[2] - facet[0];
	CrossProduct(p1, p2, normal);

	Normalize(normal);
	return normal;
}

void CopyPoints(Point3f *points, Point3f *result, int size)
{
	for (int i = 0; i <= size; ++i)
	{
		result[i].cx = points[i].cx;
		result[i].cy = points[i].cy;
		result[i].cz = points[i].cz;
	}
}

Point3f CenterOfPolygon(const Point3f *facet, int size)
{
	Point3f p(0, 0, 0);

	for (int i = 0; i < size; ++i)
	{
		p = p + facet[i];
	}

	return p/size;
}

