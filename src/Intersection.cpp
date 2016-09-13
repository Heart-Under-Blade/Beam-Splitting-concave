#include "Intersection.h"

bool inside(const Point3f &x, const Point3f &p1, const Point3f &p2, const Point3f &normal)
{
	__m128 _x = _mm_load_ps(x.point);
	__m128 _p1 = _mm_load_ps(p1.point);
	__m128 _p2 = _mm_load_ps(p2.point);
	__m128 _n = _mm_load_ps(normal.point);

	__m128 eps = _mm_set_ss(EPS_INTERSECTION);
	__m128 m_eps = _mm_set_ss(-EPS_INTERSECTION);

	__m128 a = _mm_sub_ps(_p2, _p1);
	__m128 b = _mm_sub_ps(_x, _p1);
	__m128 dir = _mm_dp_ps(_cross_product(a, b), _n, MASK_1LOW);

	return _mm_ucomigt_ss(dir, m_eps) && _mm_ucomilt_ss(dir, eps);
}

void computeIntersection(const Point3f &s, const Point3f &e,
						 const Point3f &p1, const Point3f &p2,
						 const Point3f &normal, Point3f &x)
{
	__m128 _a = _mm_load_ps(s.point);
	__m128 _b = _mm_load_ps(e.point);

	__m128 _p1 = _mm_load_ps(p1.point);
	__m128 _p2 = _mm_load_ps(p2.point);

	__m128 n = _mm_load_ps(normal.point);

	__m128 ab = _mm_sub_ps(_b, _a);
	__m128 p12 = _mm_sub_ps(_p2, _p1);

	/// normal of new plane
	__m128 n2 = _cross_product(p12, n);

	/// normalize normal
	__m128 nn2 = _mm_mul_ps(n2, _mm_rsqrt_ps(_mm_dp_ps(n2, n2, MASK_FULL)));

	/// D
	__m128 m_d_param = _mm_dp_ps(_p1, nn2, MASK_FULL);

	/// intersection vector and new plane
	__m128 dp0 = _mm_dp_ps(ab, nn2, MASK_FULL);
	__m128 dp1 = _mm_dp_ps(_a, nn2, MASK_FULL);

	__m128 add = _mm_sub_ps(dp1, m_d_param);
	__m128 t = _mm_div_ps(add, dp0);

	__m128 mul = _mm_mul_ps(t, ab);
	__m128 p = _mm_sub_ps(_a, mul);

	x.point[0] = p[0];
	x.point[1] = p[1];
	x.point[2] = p[3];
}
