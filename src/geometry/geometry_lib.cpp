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

double NormD(const Point3d &p)
{
	return	  p.x * p.x
			+ p.y * p.y
			+ p.z * p.z;
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

// REF: try to move to Point3f
Point3f CrossProduct(const Point3f &v1, const Point3f &v2)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	Point3f res;
	res.cx = _cp[0];
	res.cy = _cp[1];
	res.cz = _cp[2];

	return res;
}

// OPT:
Point3d CrossProductD(const Point3d &v1, const Point3d &v2)
{
	Point3d res;
	res.x = v1.y*v2.z - v1.z*v2.y;
	res.y = v1.z*v2.x - v1.x*v2.z;
	res.z = v1.x*v2.y - v1.y*v2.x;
	return res;
}

double Length(const Point3f &v)
{
	return sqrt(Norm(v));
}

double LengthD(const Point3d &v)
{
	return sqrt(NormD(v));
}

void Normalize(Point3f &v)
{
	double lenght = sqrt(Norm(v));
	v.cx /= lenght;
	v.cy /= lenght;
	v.cz /= lenght;
}
