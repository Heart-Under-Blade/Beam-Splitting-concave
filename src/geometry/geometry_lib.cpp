#include <math.h>
#include "Particle.h"
#include "Intersection.h"
#include "intrinsic/intrinsics.h"

float DotProduct(const Vector3f &v1, const Vector3f &v2)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _dp0 = _mm_dp_ps(_v1, _v2, MASK_FULL);
	return _dp0[0];
}

double DotProductD(const Vector3d &v1, const Vector3d &v2)
{
	return	  v1.x * v2.x
			+ v1.y * v2.y
			+ v1.z * v2.z;
}

double Norm(const Vector3f &p)
{
	return	  p.cx * p.cx
			+ p.cy * p.cy
			+ p.cz * p.cz;
}

double NormD(const Vector3d &p)
{
	return	  p.x * p.x
			+ p.y * p.y
			+ p.z * p.z;
}

void CrossProduct(const Vector3f &v1, const Vector3f &v2, Vector3f &res)
{
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _cp = _cross_product(_v1, _v2);

	res.cx = _cp[0];
	res.cy = _cp[1];
	res.cz = _cp[2];
}

Point3f IntersectVectors(const Point3f &c1, const Point3f &c2,
						 const Point3f &v1, const Point3f &v2,
						 const Point3f &normalToFacet, bool &isOk)
{
	__m128 _c1 = _mm_setr_ps(c1.cx, c1.cy, c1.cz, 0.0);
	__m128 _c2 = _mm_setr_ps(c2.cx, c2.cy, c2.cz, 0.0);
	__m128 _v1 = _mm_setr_ps(v1.cx, v1.cy, v1.cz, 0.0);
	__m128 _v2 = _mm_setr_ps(v2.cx, v2.cy, v2.cz, 0.0);
	__m128 _en = _mm_setr_ps(normalToFacet.cx,
							 normalToFacet.cy,
							 normalToFacet.cz, 0.0);

	__m128 _x = intersect_iv(_c1, _c2, _v1, _v2, _en, isOk);
	return Point3f(_x[0], _x[1], _x[2]);
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

std::ostream &operator <<(std::ostream &os, const Point3f &p)
{
	os << p.coordinates[0]
			<< " " << p.coordinates[1]
			<< " " << p.coordinates[2];
	return os;
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

double Length(const Vector3f &v)
{
	return sqrt(Norm(v));
}

double LengthD(const Vector3d &v)
{
	return sqrt(NormD(v));
}

void Normalize(Vector3f &v)
{
	double len = Length(v);
	v.cx /= len;
	v.cy /= len;
	v.cz /= len;
}

Vector3d NormalizeD(const Vector3d &v)
{
	Vector3d res;
	double lenght = sqrt(NormD(v));
	res.x = v.x / lenght;
	res.y = v.y / lenght;
	res.z = v.z / lenght;
	return res;
}

Point3f ProjectPointToPlane(const Point3f &point, const Vector3f &direction,
							const Vector3f &planeNormal)
{
	double tmp = DotProduct(point, planeNormal);
	double dp  = DotProduct(direction, planeNormal);
	tmp += planeNormal.d_param;
	tmp /= dp;
	return point - (direction * tmp);
}
