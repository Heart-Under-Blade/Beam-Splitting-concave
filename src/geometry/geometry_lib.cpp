#include <math.h>
#include "Particle.h"

#include "Intersection.h"
#include "intrinsic/intrinsics.h"

Point3f Geometry::ProjectPointToPlane(const Point3f &point,
									  const Vector3f &direction,
									  const Vector3f &planeNormal)
{
	double tmp = Point3f::DotProduct(point, planeNormal);
	double dp  = Point3f::DotProduct(direction, planeNormal);
	tmp += planeNormal.d_param;
	tmp /= dp;
	return point - (direction * tmp);
}

double Norm(const Vector3f &p)
{
	return	  p.coordinates[0] * p.coordinates[0]
			+ p.coordinates[1] * p.coordinates[1]
			+ p.coordinates[2] * p.coordinates[2];
}

Point3f IntersectVectors(const Point3f &c1, const Point3f &c2,
						 const Point3f &v1, const Point3f &v2,
						 const Point3f &normalToFacet, bool &isOk)
{
	__m128 _c1 = _mm_setr_ps(c1.coordinates[0], c1.coordinates[1], c1.coordinates[2], 0.0);
	__m128 _c2 = _mm_setr_ps(c2.coordinates[0], c2.coordinates[1], c2.coordinates[2], 0.0);
	__m128 _v1 = _mm_setr_ps(v1.coordinates[0], v1.coordinates[1], v1.coordinates[2], 0.0);
	__m128 _v2 = _mm_setr_ps(v2.coordinates[0], v2.coordinates[1], v2.coordinates[2], 0.0);
	__m128 _en = _mm_setr_ps(normalToFacet.coordinates[0],
							 normalToFacet.coordinates[1],
							 normalToFacet.coordinates[2], 0.0);

	__m128 _x = intersect_iv(_c1, _c2, _v1, _v2, _en, isOk);
	return Point3f(_x[0], _x[1], _x[2]);
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

void Geometry::DifferPolygons(const Polygon &subject, const Vector3f &subjNormal,
							  const Polygon &clip, const Vector3f &clipNormal,
							  const Vector3f &clipDir, PolygonStack &difference)
{
	__m128 _clip[MAX_VERTEX_NUM];
	bool isProjected = ProjectPolygonToPlane(clip, clipDir, subjNormal, _clip);

	if (!isProjected)
	{
		difference.Push(subject);
		return;
	}

	__m128 _clip_normal = _mm_setr_ps(clipNormal.coordinates[0], clipNormal.coordinates[1], clipNormal.coordinates[2], 0.0);

	int clipSize = clip.nVertices;
	__m128 _diff_pol[MAX_VERTEX_NUM];

	__m128 _subject[MAX_VERTEX_NUM];
	__m128 _buffer[MAX_VERTEX_NUM];

	for (int i = 0; i < subject.nVertices; ++i)
	{
		_subject[i] = _mm_load_ps(subject.vertices[i].coordinates);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.nVertices;

	__m128 _first_p, _second_p;
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (int i = 0; i < clip.nVertices; ++i)
	{
		int difSize = 0;

		__m128 _p1 = _p2;
		_p2 = _clip[i];

		__m128 *_tmp = _buff;
		_buff = _subj;
		_subj = _tmp;

		int subjSize = bufSize;
		bufSize = 0;

		_first_p = _subj[subjSize-1];
		isInFirst = is_inside_i(_first_p, _p1, _p2, _clip_normal);

		bool isIntersected;

		for (int j = 0; j < subjSize; ++j)
		{
			_second_p = _subj[j];
			isInSecond = is_inside_i(_second_p, _p1, _p2, _clip_normal);

			if (isInSecond)
			{
				if (!isInFirst)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
										   _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_diff_pol[difSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_buff[bufSize++] = _second_p;
			}
			else
			{
				if (isInFirst)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
										   _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_diff_pol[difSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_diff_pol[difSize++] = _second_p;
			}

			_first_p = _second_p;
			isInFirst = isInSecond;
		}

		if (difSize >= MIN_VERTEX_NUM)
		{
			Polygon resPolygon;
			RefineOutputPolygon(_diff_pol, difSize, resPolygon);

			if (resPolygon.nVertices >= MIN_VERTEX_NUM)
			{
				difference.Push(resPolygon);
			}
		}
	}
}

/// NOTE: вершины пучка и грани должны быть ориентированы в одном направлении
bool Geometry::IncidentBeamToFacet(Facet *facet, const Polygon &beamPol,
								   bool isInside, const Vector3f &incDir,
								   Polygon &intersection)
{
	const Vector3f &facetNormal = isInside ? facet->in_normal : facet->ex_normal;

	__m128 _output_points[MAX_VERTEX_NUM];
	// REF: перенести в случай невыпуклых частиц
	const Point3f &normal = facet->in_normal;
	const Point3f &normal1 = facetNormal;

	bool isProjected = ProjectPolygonToPlane(beamPol, incDir, normal1,
											 _output_points);
	if (!isProjected)
	{
		return false;
	}

	__m128 _normal_to_facet = _mm_setr_ps(-normal.coordinates[0], -normal.coordinates[1], -normal.coordinates[2], 0.0);
	__m128 *_output_ptr = _output_points;
	int outputSize = beamPol.nVertices;

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	__m128 _p1, _p2; // vertices of facet
	__m128 _s_point, _e_point;	// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = facet->vertices[facet->nVertices-1];
	_p2 = _mm_load_ps(p2.coordinates);

	for (int i = 0; i < facet->nVertices; ++i)
	{
		_p1 = _p2;
		p2 = facet->vertices[i];
		_p2 = _mm_load_ps(p2.coordinates);

		bufferSize = outputSize;
		outputSize = 0;

		__m128 *_temp = _output_ptr;
		_output_ptr = _buffer_ptr;
		_buffer_ptr = _temp;

		_s_point = _buffer_ptr[bufferSize-1];
		isInsideS = is_inside_i(_s_point, _p1, _p2, _normal_to_facet);

		bool isIntersected;

		for (int j = 0; j < bufferSize; ++j)
		{
			_e_point = _buffer_ptr[j];
			isInsideE = is_inside_i(_e_point, _p1, _p2, _normal_to_facet);

			if (isInsideE)
			{
				if (!isInsideS)
				{
					__m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
										   _normal_to_facet, isIntersected);
					if (isIntersected)
					{
						_output_ptr[outputSize++] = x;
					}
				}

				_output_ptr[outputSize++] = _e_point;
			}
			else if (isInsideS)
			{
				__m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
									   _normal_to_facet, isIntersected);
				if (isIntersected)
				{
					_output_ptr[outputSize++] = x;
				}
			}

			_s_point = _e_point;
			isInsideS = isInsideE;
		}
	}

	RefineOutputPolygon(_output_ptr, outputSize, intersection);
	return intersection.nVertices >= MIN_VERTEX_NUM;
}

bool Geometry::ProjectPolygonToPlane(const Polygon &polygon, const Vector3f &dir,
									 const Point3f &normal, __m128 *_projection)
{
	__m128 _normal = _mm_setr_ps(normal.coordinates[0], normal.coordinates[1], normal.coordinates[2], 0.0);
	__m128 _direction = _mm_setr_ps(dir.coordinates[0], dir.coordinates[1], dir.coordinates[2], 0.0);

	__m128 _d_param = _mm_set_ps1(normal.d_param);
	__m128 _dp0 = _mm_dp_ps(_direction, _normal, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_PROJECTION)
	{
		return false; /// beam is parallel to facet
	}

	for (int i = 0; i < polygon.nVertices; ++i)
	{
		const Point3f &p = polygon.vertices[i];
		__m128 _point = _mm_setr_ps(p.coordinates[0], p.coordinates[1], p.coordinates[2], 0.0);
		__m128 _dp1 = _mm_dp_ps(_point, _normal, MASK_FULL);
		__m128 _add = _mm_add_ps(_dp1, _d_param);
		__m128 _t = _mm_div_ps(_add, _dp0);
		__m128 _mul = _mm_mul_ps(_t, _direction);

		_projection[i] = _mm_sub_ps(_point, _mul);
	}

	return true;
}

void Geometry::RefineOutputPolygon(__m128 *_output_points, int outputSize,
								   Polygon &polygon)
{
	Point3f p;

	__m128 eps = _mm_load_ps1(&EPS_MERGE);
	__m128 sign_mask = _mm_set1_ps(-0.f);

	__m128 p0 = _output_points[outputSize-1];

	for (int i = 0; i < outputSize; ++i)
	{
		__m128 difference = _mm_sub_ps(_output_points[i], p0);
		__m128 abs = _mm_andnot_ps(sign_mask, difference);
		__m128 cmp = _mm_cmplt_ps(eps, abs);

		int res = _mm_movemask_ps(cmp) & 0b111;

		if (res != 0)
		{
			p.coordinates[0] = _output_points[i][0];
			p.coordinates[1] = _output_points[i][1];
			p.coordinates[2] = _output_points[i][2];
			polygon.vertices[polygon.nVertices++] = p;
		}

		p0 = _output_points[i];
	}
}
