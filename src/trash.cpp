/**
  ЗДЕСЬ БУДУТ ХРАНИТЬСЯ ФУНКЦИИ, КОТОРЫЕ МОГУТ ПРИГОДИТЬСЯ ПОЗЖЕ
  */

#include "Tracing.h"
#include "BeamClipper.h"
#include "clipper.hpp"
#include "macro.h"

typedef std::vector<std::vector<Point3f>> Polygons;

using namespace ClipperLib;

void FindDividePoint(const std::vector<Point3f> &polygon,
									 int i0, int i1, const Point3f &normal,
									 Point3f &x, int &nextPointIndex)
{
	int size = polygon.size();
	__m128 res;

	__m128 _a0 = _mm_setr_ps(polygon[i0].cx, polygon[i0].cy, polygon[i0].cz, 0.0);
	__m128 _a1 = _mm_setr_ps(polygon[i1].cx, polygon[i1].cy, polygon[i1].cz, 0.0);
	__m128 _n = _mm_load_ps(normal.point);

	bool isOk = false;
	int j = size-1;

	for (int i = 0; i < size; ++i)
	{
		__m128 _b0 = _mm_setr_ps(polygon[j].cx, polygon[j].cy, polygon[j].cz, 0.0);
		__m128 _b1 = _mm_setr_ps(polygon[i].cx, polygon[i].cy, polygon[i].cz, 0.0);

		if (is_inside_i(_b1, _a0, _a1, _n)
				&& !is_inside_i(_b0, _a0, _a1, _n))
		{
			res = intersect_i(_a0, _a1, _b0, _b1, _n, isOk); // OPT: написать опт. вариант (с уже готовыми векторами v1, v2)

			if (isOk)
			{
				x.point[0] = res[0];
				x.point[1] = res[1];
				x.point[2] = res[2];
				nextPointIndex = i;
				return;
			}
		}

		j = i;
	}

	LOG_ASSERT(false && "Divide point is not found");
}

void FillSubpolygon(int begin, int end,
									const std::vector<Point3f> &polygon,
									std::vector<Point3f> &subpolygon)
{
	for (int j = begin; j != end; ++j)
	{
		if (j == (int)polygon.size())
		{
			j = -1;
			continue;
		}

		subpolygon.push_back(polygon[j]);
	}
}

void DividePolygon(const std::vector<Point3f> &polygon, const Point3f &normal, Polygons &polygons)
{
	int size = polygon.size();
	int baseI = size-1;

	for (int i = 1; i < size; ++i)
	{
		Point3f v1 = polygon[i-1] - polygon[baseI];
		Point3f v2 = polygon[i] - polygon[baseI];
		Point3f res;
		CrossProduct(v1, v2, res);
		double dir = DotProduct(res, normal);

		if (dir < -EPS_INTERSECTION) // cavity is here
		{
			Point3f x;
			int next;
			FindDividePoint(polygon, baseI, i-1, normal,
							x, next);

			auto Divide = [&] (int begin, int end)
			{
				std::vector<Point3f> subpolygon;
				subpolygon.push_back(polygon[end]);
				subpolygon.push_back(x);
				FillSubpolygon(begin, end, polygon, subpolygon);
				DividePolygon(subpolygon, normal, polygons);
			};

			Divide(next, i-1);
			Divide(i-1, next);

			return;
		}

		baseI = i-1;
	}

	polygons.push_back(polygon);
}

void DivideConcavePolygon(const Point3f *polygon, int size,
										  const Point3f &normal,
										  Polygons &polygons)
{
	std::vector<Point3f> vec_polygon;

	for (int i = 0; i < size; ++i)
	{
		vec_polygon.push_back(polygon[i]);
	}

	DividePolygon(vec_polygon, normal, polygons); // recursive
}

//bool isOrderReversed(const Point3f oldNormal, const Path polygon)
//{
//	Point3f facet[MAX_VERTEX_NUM];

//	for (int i = 0; i < 3; ++i)
//	{
//		facet[i].cx = (float)polygon[i].X / MULTI_INDEX;
//		facet[i].cy = (float)polygon[i].Y / MULTI_INDEX;
//		facet[i].cz = (float)polygon[i].Z / MULTI_INDEX;
//	}

//	Point3f newNormal = NormalToFacet(facet);
//	double cosNO = DotProduct(newNormal, oldNormal);
//	return (cosNO < 0);
//}

void InversePolygonOrder(Path &polygon)
{
	IntPoint tmp;
	int size = polygon.size()-1;

	for (int vi = 0; vi <= size/2; ++vi)
	{
		tmp = polygon[vi];
		polygon[vi] = polygon[size-vi];
		polygon[size-vi] = tmp;
	}
}

