#pragma once

#include "intrinsic/intrinsics.h"

#define EPS_PROJECTION		0.00174532836589830883577820272085
const float EPS_INTERSECTION = 0.01;

bool inside(const Point3f &x, const Point3f &p1, const Point3f &p2, const Point3f &normal);

inline bool is_inside_i(__m128 x, __m128 p1, __m128 p2, __m128 normal)
{
	__m128 m_eps = _mm_set_ss(-EPS_INTERSECTION);

	__m128 p1_p2 = _mm_sub_ps(p2, p1);
	__m128 p1_x = _mm_sub_ps(x, p1);
	__m128 dir = _mm_dp_ps(_cross_product(p1_p2, p1_x), normal, MASK_1LOW);

	return _mm_ucomigt_ss(dir, m_eps);
}

inline __m128 computeIntersection_i(__m128 _a1, __m128 _a2, __m128 _b1, __m128 _b2, __m128 _normal_to_facet, bool &ok)
{
	__m128 _v_a = _mm_sub_ps(_a2, _a1);
	__m128 _v_b = _mm_sub_ps(_b2, _b1);

	// normal of new plane
	__m128 _normal_to_line = _cross_product(_v_b, _normal_to_facet);

	// normalize normal
	__m128 _normal_n = _normalize(_normal_to_line);

	// intersection vector and new plane
	__m128 _dp0 = _mm_dp_ps(_v_a, _normal_n, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_INTERSECTION)
	{
		ok = false;
		return _dp0;
	}

	__m128 _dp1 = _mm_dp_ps(_a1, _normal_n, MASK_FULL);
	__m128 m_d_param = _mm_dp_ps(_b1, _normal_n, MASK_FULL);

	__m128 _add = _mm_sub_ps(_dp1, m_d_param);
	__m128 _t = _mm_div_ps(_add, _dp0);

	__m128 _mul = _mm_mul_ps(_t, _v_a);

	ok = true;
	return _mm_sub_ps(_a1, _mul);
}

void computeIntersection(const Point3f &s, const Point3f &e,
						 const Point3f &p1, const Point3f &p2, const Point3f &normal,
						 Point3f &x);
