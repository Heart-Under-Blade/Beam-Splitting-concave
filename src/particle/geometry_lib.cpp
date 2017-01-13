#include <math.h>
#include "Particle.h"
#include "intrinsic/intrinsics.h"

float DotProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _dp0 = _mm_dp_ps(_v1, _v2, MASK_FULL);
	return _dp0[0];
}

double DotProductD(const Point3d &v1, const Point3d &v2)
{
	return	  v1.x * v2.x
			+ v1.y * v2.y
			+ v1.z * v2.z;
}

double Norm(const Point3f &p)
{
	return	  p.cx * p.cx
			+ p.cy * p.cy
			+ p.cz * p.cz;
}

void CrossProduct(const Point3f &v1, const Point3f &v2, Point3f &res)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	res.cx = _cp[0];
	res.cy = _cp[1];
	res.cz = _cp[2];
}

double Length(const Point3f &v)
{
	return sqrt(Norm(v));
}

void Normalize(Point3f &v)
{
	double lenght = sqrt(Norm(v));
	v.cx /= lenght;
	v.cy /= lenght;
	v.cz /= lenght;
}

Point3f NormalToPolygon(const Point3f *facet)
{
	Point3f normal;

	Point3f p1 = facet[1] - facet[0];
	Point3f p2 = facet[2] - facet[0];
	CrossProduct(p1, p2, normal);

	Normalize(normal);
	return normal;
}

Point3f CenterOfPolygon(const Polygon &polygon)
{
	Point3f p(0, 0, 0);

	for (int i = 0; i < polygon.size; ++i)
	{
		p = p + polygon.arr[i];
	}

	return p/polygon.size;
}


double AreaOfPolygon(const Polygon &p)
{
	double square = 0;
	const Point3f &basePoint = p.arr[0];
	Point3f p1 = p.arr[1] - basePoint;

	for (int i = 2; i < p.size; ++i)
	{
		Point3f p2 = p.arr[i] - basePoint;
		Point3f res;
		CrossProduct(p1, p2, res);
		square += sqrt(Norm(res));
		p1 = p2;
	}

	return square / 2.0;
}
