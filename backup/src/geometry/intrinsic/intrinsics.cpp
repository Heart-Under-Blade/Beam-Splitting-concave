#include "intrinsics.h"

//bool inLine_i(const Point3f &x, const Point3f &a, const Point3f &b)
//{
//	__m128 _x = _mm_load_ps(x.point);
//	__m128 _a = _mm_load_ps(a.point);
//	__m128 _b = _mm_load_ps(b.point);

//	__m128 ab = _mm_sub_ps(_a, _b);
//	__m128 ax = _mm_sub_ps(_a, _x);
//	__m128 bx = _mm_sub_ps(_b, _x);

//	__m128 sqr_len_ab = _mm_dp_ps(ab, ab, MASK_1LOW);
//	__m128 sqr_len_ax = _mm_dp_ps(ax, ax, MASK_1LOW);
//	__m128 sqr_len_bx = _mm_dp_ps(bx, bx, MASK_1LOW);

//	return (sqr_len_ax[0] + sqr_len_bx[0] < sqr_len_ab[0] + EPS_IN_LINE);
//}


//__m128 intersection_i(const Point3f &source_point, const Point3f &source_vector, float d_param, const Point3f &normal_vector)
//{
//	__m128 sp = _mm_load_ps(source_point.point);
//	__m128 sv = _mm_load_ps(source_vector.point);
//	__m128 nv = _mm_load_ps(normal_vector.point);
//	__m128 d = _mm_load_ps1(&d_param);

//	__m128 dp0 = _mm_dp_ps(sv, nv, MASK_FULL);
//	/// NOTE: проверить случай, когда нормаль и вектор перпендикулярны (dp0 == 0)

//	__m128 dp1 = _mm_dp_ps(sp, nv, MASK_FULL);

//	__m128 add = _mm_add_ps(dp1, d);
//	__m128 t = _mm_div_ps(add, dp0);

//	__m128 mul = _mm_mul_ps(t, sv);
//	return _mm_sub_ps(sp, mul);
//}

//void inFacet_i(const Point3f &x, const Facet &facet, PointPosition &pos)
//{
//	pos.facet_side_index_1 = -1;
//	pos.facet_side_index_2 = -1;
//	pos.position = 1;

//	int size = facet.size;

//	__m128 _x = _mm_load_ps(x.point);
//	__m128 _n = _mm_load_ps(facet.normal.point);

//	__m128 dir;
//	__m128 p1, p2;

//	__m128 eps = _mm_set_ss(EPS_IN_POLYGON);
//	__m128 m_eps = _mm_set_ss(-EPS_IN_POLYGON);

//	p2 = _mm_load_ps(facet.vertices[size-1].point);

//	for (int i = 0; i < size; ++i)
//	{
//		p1 = p2;
//		p2 = _mm_load_ps(facet.vertices[i].point);

//		__m128 a = _mm_sub_ps(p2, p1);
//		__m128 b = _mm_sub_ps(_x, p1);

//		dir = _mm_dp_ps(_cross_product(a, b), _n, 0x71);

//		if (_mm_ucomilt_ss(dir, m_eps))
//		{
//			pos.position = -1;
//			return;
//		}
//		else if (_mm_ucomilt_ss(dir, eps))
//		{
//			if (pos.position != 0)
//			{
//				pos.position = 0;
//				pos.facet_side_index_1 = i-1;
//			}
//			else {
//				pos.facet_side_index_2 = i-1;
//			}
//		}
//	}
//}

//__m128 intersection_line_facet_i(const Point3f &a, const Point3f &b, const Facet &facet)
//{
//	int size = facet.size;

//	__m128 p1, p2;
//	__m128 n = _mm_load_ps(facet.normal.point);

//	__m128 _a = _mm_load_ps(a.point);
//	__m128 _b = _mm_load_ps(b.point);
//	__m128 source_vector = _mm_sub_ps(_b, _a);

//	__m128 eps = _mm_set_ss(EPS_IN_LINE);

//	p2 = _mm_load_ps(facet.vertices[size-1].point);

//	for (int i = 0; i < size; ++i)
//	{
//		p1 = p2;
//		p2 = _mm_load_ps(facet.vertices[i].point);

//		__m128 p12 = _mm_sub_ps(p2, p1);

//		/// normal of new plane
//		__m128 n2 = _cross_product(p12, n);

//		/// normalize normal
//		__m128 nn2 = _mm_mul_ps(n2, _mm_rsqrt_ps(_mm_dp_ps(n2, n2, MASK_FULL)));

//		/// D
//		__m128 m_d_param = _mm_dp_ps(p1, nn2, MASK_FULL);

//		/// intersection vector and new plane
//		__m128 dp0 = _mm_dp_ps(source_vector, nn2, MASK_FULL);
//		/// NOTE: проверить случай, когда нормаль и вектор перпендикулярны (dp0 == 0)

//		__m128 dp1 = _mm_dp_ps(_a, nn2, MASK_FULL);

//		__m128 add = _mm_sub_ps(dp1, m_d_param);
//		__m128 t = _mm_div_ps(add, dp0);

//		__m128 mul = _mm_mul_ps(t, source_vector);
//		__m128 p = _mm_sub_ps(_a, mul);

//		/// is point in line (p1, p2)
//		__m128 ab = _mm_sub_ps(p1, p2);
//		__m128 ax = _mm_sub_ps(p1, p);
//		__m128 bx = _mm_sub_ps(p2, p);

//		__m128 sqr_len_ab = _mm_dp_ps(ab, ab, MASK_1LOW);
//		__m128 sqr_len_ax = _mm_dp_ps(ax, ax, MASK_1LOW);
//		__m128 sqr_len_bx = _mm_dp_ps(bx, bx, MASK_1LOW);

//		if (_mm_add_ss(sqr_len_ax, sqr_len_bx)[0] < _mm_add_ss(sqr_len_ab, eps)[0])
//		{
//			/// is point in line (a, b)
//			ab = _mm_sub_ps(_a, _b);
//			ax = _mm_sub_ps(_a, p);
//			bx = _mm_sub_ps(_b, p);

//			sqr_len_ab = _mm_dp_ps(ab, ab, MASK_1LOW);
//			sqr_len_ax = _mm_dp_ps(ax, ax, MASK_1LOW);
//			sqr_len_bx = _mm_dp_ps(bx, bx, MASK_1LOW);

//			if (_mm_add_ss(sqr_len_ax, sqr_len_bx)[0] < _mm_add_ss(sqr_len_ab, eps)[0])
//			{
//				return p;
//			}
//		}
//	}

//	return n;
//}

