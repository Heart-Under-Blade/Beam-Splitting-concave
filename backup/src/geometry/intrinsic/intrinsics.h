#pragma once

#include <xmmintrin.h>
#include <smmintrin.h>
#include <vector>

//#define EPS_IN_LINE 0.001
#define EPS_IN_LINE 1
#define EPS_IN_POLYGON 0.001

#define MASK_FULL 0x7F /// mask for filling __m128 - all 4 float parts
#define MASK_1LOW 0x71 /// mask for filling __m128 - only first (low) float part

#define _cross_product(a, b) _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),\
													_mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),\
										_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),\
													_mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1))))

#define _normalize(n) _mm_mul_ps(n, _mm_rsqrt_ps(_mm_dp_ps(n, n, MASK_FULL)));

/**
 * @brief The PointPosition struct
 * Position of point on the facet (uses in 'inFacet' function)
 */
struct PointPosition {
	int position; /// 1 - inside, 0 - on the side or vertex, -1 - outside
	int facet_side_index_1; /// index of side that point belongs
	int facet_side_index_2; /// index of second side (point belongs to vertex between the sides)
};

/**
 * @brief inLine_i Is point belongs to line
 * @param x Point
 * @param a First point of line
 * @param b Second point of line
 * @return True if belongs
 */
//bool inLine_i(const Point3f &x, const Point3f &a, const Point3f &b);

/** Dot product */
//__m128 res = _mm_dp_ps(_mm_load_ps(p1.point), _mm_load_ps(p2.point), 0x71);

//inline __m128 crossProduct_i(__m128 a, __m128 b)
//{
//	return _mm_sub_ps(_mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,0,2,1)),
//								 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,1,0,2))),
//					  _mm_mul_ps(_mm_shuffle_ps(a, a, _MM_SHUFFLE(3,1,0,2)),
//								 _mm_shuffle_ps(b, b, _MM_SHUFFLE(3,0,2,1))));
//}

/**
 * @brief intersection_i Intersection of vector and plane
 * @param source_point Start point of source vector
 * @param source_vector Source vector
 * @param d_param D-param of plane
 * @param normal_vector Normal of plane
 * @return Intersection point
 */
//__m128 intersection_i(const Point3f &source_point, const Point3f &source_vector,
//					  float d_param, const Point3f &normal_vector);

/**
 * @brief inFacet_i Position of point on the facet
 * @param x Point
 * @param facet Facet
 * @param pos Position of point (result)
 */
//void inFacet_i(const Point3f &x, const Facet &facet, PointPosition &pos);

//__m128 intersection_line_facet_i(const Point3f &a, const Point3f &b, const Facet &facet);
