#include "BeamClipper.h"
#include "macro.h"
#include "geometry_lib.h"
#include <cmath>

#define SWAP_COORD_THRESHOLD 0.5

using namespace ClipperLib;

BeamClipper::BeamClipper()
{
	m_clipper.ZFillFunction(FindZCoord);
}

void FindZCoord(IntPoint &a1, IntPoint &a2, IntPoint &, IntPoint &, IntPoint &point)
{
	IntPoint *top, *bot;

	if (a2.Z > a1.Z)
	{
		top = &a2;
		bot = &a1;
	}
	else
	{
		top = &a1;
		bot = &a2;
	}

	double x, y;

	x = top->X - bot->X;
	y = top->Y - bot->Y;
	double normVec = x*x + y*y;

	x = point.X - bot->X;
	y = point.Y - bot->Y;
	double normBotVec = x*x + y*y;

	point.Z = (top->Z - bot->Z) * sqrt(normBotVec) / sqrt(normVec) + bot->Z;
}

Axis BeamClipper::GetSwapAxis(const Point3f &normal) const
{
	if (fabs(normal.cx) > SWAP_COORD_THRESHOLD)
	{
		return Axis::aX;
	}
	else if (fabs(normal.cy) > SWAP_COORD_THRESHOLD)
	{
		return Axis::aY;
	}
	else // no need to swap
	{
		return Axis::aZ;
	}
}

void BeamClipper::SwapCoords(Axis ax1, Axis ax2, ClipperLib::Paths &origin) const
{
	if (ax1 == ax2)
	{
		return;
	}

	cInt *oldP = nullptr, *newP = nullptr;

	for (Path &path : origin)
	{
		for (IntPoint &p : path)
		{
			switch (ax1) {
			case Axis::aX: oldP = &p.X;
				break;
			case Axis::aY: oldP = &p.Y;
				break;
			case Axis::aZ: oldP = &p.Z;
				break;
			}

			switch (ax2) {
			case Axis::aX: newP = &p.X;
				break;
			case Axis::aY: newP = &p.Y;
				break;
			case Axis::aZ: newP = &p.Z;
				break;
			}

			cInt oldP1 = *oldP;
			cInt newP1 = *newP;
			*oldP = newP1;
			*newP = oldP1;
		}
	}
}

void BeamClipper::PolygonToPath(const Polygon &pol, Paths &path) const
{
	for (int i = 0; i < pol.size; ++i)
	{
		const Point3f &p = pol.arr[i];

		path[0] << IntPoint((cInt)(p.cx * MULTI_INDEX),
							(cInt)(p.cy * MULTI_INDEX),
							(cInt)(p.cz * MULTI_INDEX));
	}
}

double BeamClipper::AreaOfConcavePolygon(const Polygon &beam, const Point3f &normal) const
{
	double area = 0;

	Paths polygon(1);
	PolygonToPath(beam, polygon);

	Point3f normal1 = normal;

	Point3f n_normal = normal1;
	Normalize(n_normal);

	if (fabs(n_normal.cx) > SWAP_COORD_THRESHOLD)
	{
		SwapCoords(Axis::aX, Axis::aZ, polygon);
		float tmp = normal1.cx;
		normal1.cx = normal1.cz;
		normal1.cz = tmp;
	}
	else if (fabs(n_normal.cy) > SWAP_COORD_THRESHOLD)
	{
		SwapCoords(Axis::aY, Axis::aZ, polygon);
		float tmp = normal1.cy;
		normal1.cy = normal1.cz;
		normal1.cz = tmp;
	}

	area = ClipperLib::Area(polygon[0]);
	area /= (cInt)MULTI_INDEX * MULTI_INDEX;

	Point3f z(0, 0, 1);
	double cosNormZ = DotProduct(normal1, z);

	if (cosNormZ > 0)
	{
		cosNormZ = DotProduct(normal1, -z);
	}

	area = (area*sqrt(Norm(normal1)))/cosNormZ;
//	Polygons polygons;
//	DivideConcavePolygon(beam.shape, beam.shapeSize, normal,
//						 polygons);

//	for (unsigned int i = 0; i < polygons.size(); ++i)
//	{
//		square += SquareOfPolygon(polygons[i]);
//	}

	if (area < 0)
	{
		area = -area;
	}

	return area;
}

//void BeamClipper::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
//{
//	Point3f &facetNormal = m_facets[beam.facetId].ex_normal;

//	IntArray facetIds;
//	SelectVisibleFacets(beam, facetIds);

//	Paths originPath(1);
//	m_clipper.PolygonToPath(beam.polygon, originPath);

//	// cut facet projections out of beam one by one
//	for (int i = 0; i < facetIds.size; ++i)
//	{
//		int id = facetIds.arr[i];

//		Paths clippedPath(1);
//		m_clipper.CutBeamByPolygon(originPath, m_facets[id].polygon, beam.direction,
//								   facetNormal, clippedPath);

//		if (clippedPath.empty()) // beam incedents on facet totaly
//		{
//			originPath.clear();
//			break;
//		}
//		else
//		{
//			originPath = clippedPath;
//		}
//	}

//	Beam tmp = beam;

//	for (const Path &p : originPath)
//	{
//		m_clipper.PathToPolygon(p, tmp.polygon);
//		scatteredBeams.push_back(tmp);
//	}
//}

void BeamClipper::ProjectPointToFacet(const Point3d &point, const Point3d &direction,
										 const Point3d &facetNormal, Point3d &projection)
{
	double tmp = DotProductD(point, facetNormal);
	tmp = tmp + facetNormal.d;
	double dp = DotProductD(direction, facetNormal);
	tmp = tmp/dp;
	projection = point - (direction * tmp);
}

void BeamClipper::ProjectFacetToFacet(const Polygon &a_facet, const Point3f &a_dir,
									  const Point3f &b_normal, Path &projection)
{
	for (int i = 0; i < a_facet.size; ++i)
	{
		Point3d p;
		ProjectPointToFacet(Point3d(a_facet.arr[i]), Point3d(a_dir),
							Point3d(b_normal), p);

		projection << IntPoint((cInt)(p.x * MULTI_INDEX),
							   (cInt)(p.y * MULTI_INDEX),
							   (cInt)(p.z * MULTI_INDEX));
	}
}

void BeamClipper::CutBeamByPolygon(Paths &beamPol, const Polygon &polygon,
								   const Point3f &direction,
								   const Point3f &polNormal, Paths &result)
{
	Axis axis = GetSwapAxis(polNormal);
	SwapCoords(axis, Axis::aZ, beamPol);

	// проецируем грань на начальный пучок
	Paths clip(1);
	ProjectFacetToFacet(polygon, direction, polNormal, clip[0]);
	SwapCoords(axis, Axis::aZ, clip);

	Difference(beamPol, clip, result);

	if (!result.empty())
	{
		HandleResultPaths(axis, result);
	}
}

void BeamClipper::Difference(const Paths &subject, const Paths &clip, Paths &difference)
{
	m_clipper.AddPaths(subject, ptSubject, true);
	m_clipper.AddPaths(clip, ptClip, true);
	m_clipper.Execute(ctDifference, difference);
	m_clipper.Clear();
}

// NOTE: функция расчитана не более чем на 2 полигона
void BeamClipper::RemoveHole(Paths &result) const
{
	Path &pol1 = result.front();
	Path &pol2 = result.back();

	bool or1 = ClipperLib::Orientation(pol1);
	bool or2 = ClipperLib::Orientation(pol2);

	if (or1 == or2)
	{
		return;
	}

	unsigned first = 0, second = 0;
	cInt len;
	cInt min_len = LONG_MAX;
	IntPoint diff;

	for (unsigned i = 0; i < pol1.size(); ++i)
	{
		for (unsigned j = 0; j < pol2.size(); ++j)
		{
			diff.X = pol1.at(i).X - pol2.at(j).X;
			diff.Y = pol1.at(i).Y - pol2.at(j).Y;
			diff.Z = pol1.at(i).Z - pol2.at(j).Z;

			len = diff.X*diff.X + diff.Y*diff.Y + diff.Z*diff.Z;
			len = (cInt)std::round(sqrt((double)len));

			if (len < min_len)
			{
				first = i;
				second = j;
				min_len = len;
			}
		}
	}

	Paths newResult(1);

	for (unsigned i = 0; i < pol1.size(); ++i)
	{
		newResult[0] << pol1.at(i);

		if (i == first)
		{
			unsigned j = second;

			do
			{
				if (j == pol2.size())
				{
					j = 0;
				}
				else
				{
					newResult[0] << pol2.at(j);
					++j;
				}
			}
			while (j != second);

			newResult[0] << pol2.at(second);
			newResult[0] << pol1.at(first);
		}
	}

	result = newResult;
}

void BeamClipper::RemoveEmptyPaths(Paths &result) const
{
	Paths buff = result;
	result.clear();

	for (unsigned i = 0; i < buff.size(); ++i)
	{
		if (buff.at(i).size() >= CLIP_RESULT_SINGLE)
		{
			result.push_back(buff.at(i));
		}
	}
}

// TODO: delete this class totally

void BeamClipper::HandleResultPaths(Axis axis, Paths &result) const
{
	SwapCoords(Axis::aZ, axis, result); // меняем координаты обратно
	ClipperLib::CleanPolygons(result, EPS_MULTI);
	RemoveEmptyPaths(result);

//	LOG_ASSERT(result.size() < 3);

	if (result.size() == 2)
	{
		RemoveHole(result);
	}
}

void BeamClipper::PathToPolygon(const Path &path, Polygon &polygon) const
{
	int vertexNum = path.size();
	polygon.size = vertexNum;

	for (const IntPoint &p : path)
	{
		Point3f tmp((double)p.X / MULTI_INDEX,
					(double)p.Y / MULTI_INDEX,
					(double)p.Z / MULTI_INDEX);

		polygon.arr[--vertexNum] = tmp;
	}

//	LOG_ASSERT(polygon.size > 0 && polygon.size <= MAX_VERTEX_NUM);
}
