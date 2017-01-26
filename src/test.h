/**
  File for various tests
*/

#pragma once

#include "global.h"
#include <fstream>
#include <iostream>

#include "Hexagonal.h"
#include "ConcaveHexagonal.h"
#include "Intersection.h"
#include "Tracing.h"
void SetOutputPolygon(__m128 *_output_points, int outputSize,
							   Polygon &polygon)
{
	Point3f p;

	__m128 eps = _mm_load_ps1(&EPS_INTERSECTION);
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
			p.cx = _output_points[i][0];
			p.cy = _output_points[i][1];
			p.cz = _output_points[i][2];
			polygon.arr[polygon.size++] = p;
		}

		p0 = _output_points[i];
	}
}
bool ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
								  const Point3f &normal, __m128 *_projection)
{
	__m128 _normal = _mm_setr_ps(normal.cx, normal.cy, normal.cz, 0.0);
	__m128 _direction = _mm_setr_ps(dir.cx, dir.cy, dir.cz, 0.0);

	__m128 _d_param = _mm_set_ps1(normal.d_param);
	__m128 _dp0 = _mm_dp_ps(_direction, _normal, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_PROJECTION)
	{
		return false; /// beam is parallel to facet
	}

	for (int i = 0; i < polygon.size; ++i)
	{
		const Point3f &p = polygon.arr[i];
		__m128 _point = _mm_setr_ps(p.cx, p.cy, p.cz, 0.0);
		__m128 _dp1 = _mm_dp_ps(_point, _normal, MASK_FULL);
		__m128 _add = _mm_add_ps(_dp1, _d_param);
		__m128 _t = _mm_div_ps(_add, _dp0);
		__m128 _mul = _mm_mul_ps(_t, _direction);

		_projection[i] = _mm_sub_ps(_point, _mul);
	}

	return true;
}

void Differ(const Polygon &clip, const Point3f &normal,
						 const Polygon &subject, const Point3f &subjectDir,
						 Polygon *difference, int &resultSize)
{
	__m128 _clip[MAX_VERTEX_NUM];
	bool isProjected = ProjectToFacetPlane(clip, subjectDir, normal, _clip);

	if (!isProjected)
	{
		difference[resultSize++] = subject;
		return;
	}

	__m128 _clip_normal = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
	int clipSize = clip.size;
	__m128 _res_pol[MAX_VERTEX_NUM]; // OPT: заменить на Polygon

	__m128 _subject[MAX_VERTEX_NUM];
	__m128 _buffer[MAX_VERTEX_NUM];

	for (int i = 0; i < subject.size; ++i)
	{
		_subject[i] = _mm_load_ps(subject.arr[i].point);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.size;

	__m128 _first_p, _second_p; // points of projection
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (int i = 0; i < clip.size; ++i)
	{
		int resSize = 0;

		__m128 _p1 = _p2;
		_p2 = _clip[i];

		__m128 *_tmp = _buff;
		_buff = _subj;
		_subj = _tmp;

		int subjSize = bufSize;
		bufSize = 0;

		_first_p = _subj[subjSize-1];
		isInSecond = is_inside_i(_first_p, _p1, _p2, _clip_normal);

		bool isIntersected;

		for (int j = 0; j < subjSize; ++j)
		{
			_second_p = _subj[j];
			isInFirst = is_inside_i(_second_p, _p1, _p2, _clip_normal);

			if (isInFirst)
			{
				if (!isInSecond)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
													 _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _p1, _second_p))
					{
						_res_pol[resSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_buff[bufSize++] = _second_p;
			}
			else
			{
				if (isInSecond)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
													 _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_res_pol[resSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_res_pol[resSize++] = _second_p;
			}

			_first_p = _second_p;
			isInSecond = isInFirst;
		}

		if (resSize >= MIN_VERTEX_NUM)
		{
			Polygon resPolygon;
			SetOutputPolygon(_res_pol, resSize, resPolygon);

			if (resPolygon.size >= MIN_VERTEX_NUM)
			{
				difference[resultSize++] = resPolygon;
			}
		}
	}
}

void testDiff()
{
	Polygon a;
	a.arr[0] = Point3f(0, 0, 0);
	a.arr[1] = Point3f(0, 4, 0);
	a.arr[2] = Point3f(4, 4, 0);
	a.arr[3] = Point3f(4, 0, 0);
	a.size = 4;

	Point3f n(0, 0, 1);

	Polygon s;
	s.arr[0] = Point3f(2, 2, 1);
	s.arr[1] = Point3f(2, 6, 1);
	s.arr[2] = Point3f(6, 6, 1);
	s.arr[3] = Point3f(6, 2, 1);
	s.size = 4;

	Polygon r[32];
	int size = 0;

	Differ(s, n, a, -n, r, size);

	int fff = 0;
}

void outputParticle(const Particle &particle)
{
	for (int i = 0; i < particle.facetNum; ++i)
	{
		std::cout << i << ": ";
		for (int j = 0; j < particle.facets[i].polygon.size; ++j)
		{
			Point3f p = particle.facets[i].polygon.arr[j];
			std::cout << "("
					  << p.point[0] << ", "
					  << p.point[1] << ", "
					  << p.point[2]
					  << "), ";
		}

		std::cout << std::endl;
	}

	std::cout << std::endl << "Normals" << std::endl << std::endl;

	for (int i = 0; i < particle.facetNum; ++i)
	{
		std::cout << i << ": ";
		std::cout << "("
				  << particle.facets[i].in_normal.cx << ", "
				  << particle.facets[i].in_normal.cy << ", "
				  << particle.facets[i].in_normal.cz
				  << "), ";
		std::cout << std::endl;
	}
}

void toFile(const Particle &particle)
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < particle.facetNum; ++i)
	{
		for (int j = 0; j < particle.facets[i].polygon.size; ++j)
		{
			Point3f p = particle.facets[i].polygon.arr[j];
			M << p.point[0] << ' '
							<< p.point[1] << ' '
							<< p.point[2];
			M << std::endl;
		}

		M << std::endl << std::endl;;
	}

	M.close();
}

void testHexagonBuilding()
{
	Particle hex = Hexagonal(20, 100, 1.31);
	toFile(hex);
	outputParticle(hex);
}

void testConcaveHexagonRot()
{
	Particle *hex = new ConcaveHexagonal(20, 100, 1.31, 10);
	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
	hex->Rotate(beta, gamma, 0);
	toFile(*hex);
	outputParticle(*hex);
}

void testHexagonRotate()
{
	Particle *hex = new Hexagonal(20, 100, 1.31);
//	double beta = 19*M_PI/180;
//	double gamma = 90*M_PI/180;
	double beta = ((19 + 0.5)*M_PI)/(2.0*100);
	double gamma = ((36 + 0.5)*M_PI)/(3.0*101);
	hex->Rotate(beta, gamma, 0);

	toFile(*hex);
	outputParticle(*hex);
}


