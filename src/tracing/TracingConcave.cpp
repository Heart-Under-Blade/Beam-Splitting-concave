#include "TracingConcave.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>

#ifdef _DEBUG // DEB
#include <iostream>
#endif

#define EPS_ORTO_FACET 0.0001 //TODO: подобрать норм значение

using namespace ClipperLib;

#ifdef _TRACK_ALLOW
std::ofstream trackMapFile("tracks.dat", std::ios::out);
#endif

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	Point3f point = m_initialBeam.direction * m_particle->GetHalfHeight()*2;
	m_initialBeam.direction.d_param = DotProduct(point, m_initialBeam.direction);
	m_initialBeam.location = Location::Outside;

	m_isArea = false;
}

void TracingConcave::SplitBeamByParticle(std::vector<Beam> &outBeams)
{
	m_treeSize = 0;

	TraceFirstBeam();
	TraceSecondaryBeams(outBeams);

#ifdef _DEBUG // DEB
//	double rrr = 0;
//	for (const Beam &b : outBeams)
//	{
//		if (b.track.size() == 1 &&
//				b.track[0] == 1)
//			rrr += AreaOfBeam(b);
//	}
//	int ggg = 0;
#endif
}

void TracingConcave::TraceFirstBeam()
{
#ifdef _TRACK_OUTPUT
	trackMapFile << "0 lvl: ";
#endif
	m_lightSurfaceArea = 0;

	IntArray orderedFacets;
	SelectFirstBeamVisibleFacets(orderedFacets);

	for (int i = 0; i < orderedFacets.size; ++i)
	{
		bool hasIntersection;

		Polygon resFacets[MAX_POLYGON_NUM];
		int resSize;
		IntersectWithFacet(orderedFacets, i, resFacets, resSize, hasIntersection);

		if (hasIntersection)
		{
			Beam inBeam, outBeam;

			// set optical params of beam
			int facetId = orderedFacets.arr[i];
			SetFirstBeamOpticalParams(facetId, inBeam, outBeam);

#ifdef _DEBUG // DEB
			double eee = 0;
#endif
			for (int j = 0; j < resSize; ++j)
			{
#ifdef _DEBUG // DEB
				eee += AreaOfPolygon(resFacets[j]);
#endif
				// set geometry of beam
				 inBeam.SetPolygon(resFacets[j]);
				outBeam.SetPolygon(resFacets[j]);

				PushBeamToTree( inBeam, facetId, 0, Location::Inside );
				PushBeamToTree(outBeam, facetId, 0, Location::Outside);

				if (m_isArea)
				{
					CalcLigthSurfaceArea(facetId, outBeam);
				}
			}
		}
	}
}

void TracingConcave::SelectFirstBeamVisibleFacets(IntArray &facetIds)
{
	FindVisibleFacets_initial(facetIds);
	SortFacets(m_initialBeam.direction, facetIds);
}

void TracingConcave::IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
											Polygon *resFacets, int &resSize,
											bool &hasIntersection)
{
	resSize = 0;
	hasIntersection = true;
	int facetId = facetIds.arr[prevFacetNum];

	// OPT: раскомментить
	if (prevFacetNum == 0 /*|| m_particle->IsUnshadowedExternal(facetId)*/) // this facet is obviously not shadowed
	{
		Polygon beamPolygon;
		SetPolygonByFacet(facetId, beamPolygon);
		resFacets[resSize++] = beamPolygon;
	}
	else // facet is probably shadowed by others
	{
		CutShadowsFromFacet2(facetId, facetIds, prevFacetNum, m_initialBeam,
							 resFacets, resSize);
		if (resSize == 0)
		{
			hasIntersection = false;
		}
	}
}

void TracingConcave::SelectVisibleFacets(const Beam &beam, IntArray &facetIds)
{
	FindVisibleFacets(beam, facetIds);

	Point3f dir = beam.direction;
	dir.d_param = m_facets[beam.facetId].in_normal.d_param;
	SortFacets(dir, facetIds);
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
{
	Point3f &facetNormal = m_facets[beam.facetId].ex_normal;

	IntArray facetIds;
	SelectVisibleFacets(beam, facetIds);

	Paths originPath(1);
	m_clipper.PolygonToPath(beam.polygon, originPath);

	// cut facet projections out of beam one by one
	for (int i = 0; i < facetIds.size; ++i)
	{
		int id = facetIds.arr[i];

		Paths clippedPath(1);
		CutBeamByFacet(originPath, id, beam.direction, facetNormal, clippedPath);

		if (clippedPath.empty()) // beam incedents on facet totaly
		{
			originPath.clear();
			break;
		}
		else
		{
			originPath = clippedPath;
		}
	}

	Beam tmp = beam;

	for (const Path &p : originPath)
	{
		m_clipper.PathToPolygon(p, tmp.polygon);
		scatteredBeams.push_back(tmp);
	}
}

//void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
//{
//	Point3f &facetNormal = m_facets[beam.facetId].ex_normal;

//	IntArray facetIds;
//	SelectVisibleFacets(beam, facetIds);

//	// cut facet projections out of beam one by one
//	for (int i = 0; i < facetIds.size; ++i)
//	{
//		int id = facetIds.arr[i];

//		Polygon resultBeams[MAX_VERTEX_NUM];
//		int resultSize;

//		Difference(beam.polygon, facetNormal, m_facets[id].polygon, -beam.direction,
//				   resultBeams, resultSize);

//		if (resultSize == 0) // beam incedents on facet totaly
//		{
//			beam.polygon.size = 0;
//			break;
//		}
//		else
//		{
//			beams = clippedPath;
//		}
//	}

//	Beam tmp = beam;

//	for (const Path &p : originPath)
//	{
//		m_clipper.PathToPolygon(p, tmp.polygon);
//		scatteredBeams.push_back(tmp);
//	}
//}

void TracingConcave::PushBeamToTree(Beam &beam, int facetId, int level,
									Location location)
{
	assert(m_treeSize < MAX_BEAM_REFL_NUM);

#ifdef _TRACK_ALLOW
	AddToTrack(beam, facetId);
#ifdef _TRACK_OUTPUT
	PrintTrack(beam, facetId);
#endif
#endif
	beam.facetId = facetId;
	beam.level = level;
	beam.location = location;
	m_beamTree[m_treeSize] = beam;
	++m_treeSize;
}

#ifdef _TRACK_ALLOW
void TracingConcave::AddToTrack(Beam &beam, int facetId)
{
	std::vector<int> &tr = beam.track;
	int size = tr.size();

	if (size == 0 || (size > 0 && facetId != tr.at(size-1)))
		tr.push_back(facetId);
}

void TracingConcave::PrintTrack(const Beam &beam, int facetId)
{
	trackMapFile << "(" << beam.track.at(0);

	for (int i = 1; i < beam.track.size(); ++i)
	{
		trackMapFile << ", ";
		trackMapFile << beam.track.at(i);
	}

	if (facetId != beam.track.back())
	{
		trackMapFile << ", " << facetId;
	}

	trackMapFile << "), ";
	trackMapFile.flush();
}
#endif

void TracingConcave::CutIncidentBeam2(int facetId, Beam &beam, bool &isDivided)
{
	isDivided = false;

	// OPT: раскомментить
//	if (beam.location == Location::Inside
//			&& !m_particle->IsShadowedInternal(beam.facetId))
//	{
//		return;
//	}

	const Facet &beamFacet = m_facets[beam.facetId];
//	const Point3f &facetNormal = facet.in_normal;
	const Point3f &beamNormal = beamFacet.normal[(int)beam.location];

	const Point3f &facetNormal = (beam.location == Location::Outside) ? -beamNormal
																	  : beamNormal;

	Polygon resultBeams[MAX_VERTEX_NUM];
	int resultSize = 0;
	Difference(m_facets[facetId].polygon, facetNormal, beam.polygon, beamNormal,
			   -beam.direction, resultBeams, resultSize);

	if (resultSize == 0) // beam is totaly swallowed by facet
	{
		beam.polygon.size = 0;
	}
	else if (resultSize == CLIP_RESULT_SINGLE)
	{
		beam.polygon = resultBeams[0];
	}
	else // beam had divided by facet
	{
		Beam tmp = beam;

		for (int i = 0; i < resultSize; ++i)
		{
			tmp.polygon = resultBeams[i];
			PushBeamToTree(tmp, tmp.facetId, tmp.level, tmp.location);
		}

		isDivided = true;
		beam.polygon.size = 0;
	}
}

bool TracingConcave::HasExternalBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Outside
			&& incidentBeam.polygon.size != 0);
}

void TracingConcave::PushBeamsToTree(int level, int facetId, bool hasOutBeam,
									 Beam &inBeam, Beam &outBeam)
{
	if (hasOutBeam)
	{
		PushBeamToTree(outBeam, facetId, level, Location::Outside);
	}

#ifdef _TRACK_ALLOW
	inBeam.track = inBeam.track;
#ifdef _TRACK_OUTPUT
	trackMapFile << "[in] ";
#endif
#endif
	PushBeamToTree(inBeam, facetId, level, Location::Inside);
}

void TracingConcave::TraceSecondaryBeams(std::vector<Beam> &scaterredBeams)
{
	while (m_treeSize != 0)
	{
		Beam incidentBeam = m_beamTree[--m_treeSize];
		const Location loc = incidentBeam.location;

		if (isTerminalBeam(incidentBeam))
		{
			if (loc == Location::Outside)
			{
				CatchExternalBeam(incidentBeam, scaterredBeams);
			}

			continue;
		}

#ifdef _TRACK_OUTPUT
		trackMapFile << "\n" << incidentBeam.level << " lvl: ";
		trackMapFile.flush();
#endif
		IntArray facetIds;
		SelectVisibleFacets(incidentBeam, facetIds);

		for (int i = 0; i < facetIds.size; ++i)
		{
			int facetId = facetIds.arr[i];
			Polygon intersected;
			bool hasIntersection = Intersect(facetId, incidentBeam, intersected);

			if (hasIntersection)
			{
				Beam inBeam, outBeam;
				inBeam.SetPolygon(intersected);
				outBeam.SetPolygon(intersected);

				bool hasOutBeam;
				SetOpticalBeamParams(facetId, incidentBeam,
									 inBeam, outBeam, hasOutBeam);

				PushBeamsToTree(incidentBeam.level+1, facetId, hasOutBeam,
								inBeam, outBeam);

				bool isDivided;
				CutIncidentBeam2(facetId, incidentBeam, isDivided);

				if (isDivided)
				{
					break;
				}
			}
		}

		if (HasExternalBeam(incidentBeam))
		{	// посылаем обрезанный всеми гранями внешний пучок на сферу
			scaterredBeams.push_back(incidentBeam);
		}
	}
}

void TracingConcave::SetOpticalBeamParams(int facetId, Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam,
										  bool &hasOutBeam)
{
#ifdef _TRACK_ALLOW
	outBeam.track = incidentBeam.track;
	inBeam.track = incidentBeam.track;
#endif
	const Point3f &incidentDir = incidentBeam.direction;
	const Point3f &normal = m_facets[facetId].ex_normal;

	hasOutBeam = true;

	double cosIN = DotProduct(incidentDir, normal);

	if (cosIN >= EPS_COS_00) /// normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);
	}
	else // slopping incidence
	{
		if (incidentBeam.location == Location::Inside)
		{
			Beam incBeam = incidentBeam;
			SetSloppingIncidenceBeamParams(cosIN, normal, incBeam,
										   inBeam, outBeam, hasOutBeam);
		}
		else // beam is external
		{
			inBeam.JMatrix = incidentBeam.JMatrix;
			double cosI = DotProduct(-normal, incidentDir);

			SetSloppingBeamParams_initial(incidentDir, cosI, facetId,
										  inBeam, outBeam);
			if (m_isOpticalPath)
			{
				CalcOpticalPathInternal(cosIN, incidentBeam, inBeam, outBeam);
			}
		}
	}
}

void TracingConcave::FindVisibleFacets_initial(IntArray &facetIds)
{
	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIN = DotProduct(m_initialBeam.direction, m_facets[i].in_normal);

		if (cosIN > EPS_COS_90) // beam incidents to this facet
		{
			facetIds.arr[facetIds.size++] = i;
		}
	}
}

void TracingConcave::FindVisibleFacets(const Beam &beam, IntArray &facetIds)
{
	int type = !((bool)beam.location);
	const Point3f &invNormal = -m_facets[beam.facetId].normal[type];

	Point3f cob = CenterOfPolygon(beam.polygon);

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIN = DotProduct(beam.direction, m_facets[i].normal[type]);

		if (cosIN >= EPS_COS_90) // beam incidents to this facet
		{
			const Point3f &cof = m_particle->centers[i];
			Point3f vectorToFacet = cof - cob;
			double cosFacets = DotProduct(invNormal, vectorToFacet);

			if (cosFacets >= EPS_ORTO_FACET) // facet is in front of begin of beam
			{
				facetIds.arr[facetIds.size++] = i;
			}
		}
	}
}

void TracingConcave::CutShadowsFromFacet2(int facetId, const IntArray &facetIds,
										  int prevFacetNum, const Beam &beam,
										  Polygon *resFacets, int &resSize)
{
	const Point3f &facetNormal = m_facets[facetId].normal[(int)beam.location];
//	const Point3f &facetNormal = m_facets[facetId].in_normal;
	resFacets[resSize++] = m_facets[facetId].polygon;

	for (int i = 0; i < prevFacetNum; ++i)
	{
		int id = facetIds.arr[i];

		Polygon diffFacets[MAX_POLYGON_NUM];
		int diffSize = 0;

		while (resSize != 0)
		{
			const Polygon &clip = m_facets[id].polygon;
			const Polygon &subj = resFacets[--resSize];
//			const Point3f &clipNormal = m_facets[id].normal[(int)beam.location];
			const Point3f &clipNormal = facetNormal;

			Difference(clip, clipNormal, subj, facetNormal, -beam.direction,
					   diffFacets, diffSize);
		}

		if (diffSize != 0)
		{
			for (int i = 0; i < diffSize; ++i)
			{
				resFacets[resSize++] = diffFacets[i];
			}
		}
		else // beam is totaly swallowed by facet
		{
			break;
		}
	}
}

void TracingConcave::ProjectPointToFacet(const Point3d &point, const Point3d &direction,
										 const Point3d &facetNormal, Point3d &projection)
{
	double tmp = DotProductD(point, facetNormal);
	tmp = tmp + facetNormal.d;
	double dp = DotProductD(direction, facetNormal);
	tmp = tmp/dp;
	projection = point - (direction * tmp);
}

void TracingConcave::ProjectPointToFacet(const Point3f &point, const Point3f &direction,
										 const Point3f &facetNormal, Point3f &projection)
{
	double t = DotProduct(point, facetNormal);
	t = t + facetNormal.d_param;
	double dp = DotProduct(direction, facetNormal);
	t = t/dp;
	projection = point - (direction * t);
}

/// OPT: поменять все int и пр. параметры функций на ссылочные

void TracingConcave::ProjectFacetToFacet(const Polygon &a_facet, const Point3f &a_dir,
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

/** TODO: придумать более надёжную сортировку по близости
 * (как вариант определять, что одна грань затеняют другую по мин. и макс.
 * удалённым вершинам, типа: "//" )*/
void TracingConcave::SortFacets(const Point3f &beamDir, IntArray &facetIds)
{
	float distances[MAX_VERTEX_NUM];

	for (int i = 0; i < facetIds.size; ++i)
	{
		const int &id = facetIds.arr[i];
		distances[i] = CalcMinDistanceToFacet(m_facets[id].polygon, beamDir);
	}

	int left = 0;
	int rigth = facetIds.size - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

	while (true)
	{
		float base = distances[(left + rigth)/2];

		int i = left;
		int j = rigth;

		while (i <= j)
		{
			while (distances[i] < base)
			{
				++i;
			}

			while (distances[j] > base)
			{
				--j;
			}

			if (i <= j)	// exchange elems
			{
				float temp_d = distances[i];
				distances[i] = distances[j];
				distances[j] = temp_d;

				int temp_v = facetIds.arr[i];
				facetIds.arr[i] = facetIds.arr[j];
				facetIds.arr[j] = temp_v;

				++i;
				--j;
			}
		}

		if (i < rigth)
		{
			stack[size++] = i;
			stack[size++] = rigth;
		}

		if (left < j)
		{
			stack[size++] = left;
			stack[size++] = j;
		}

		if (size == 0)
		{
			break;
		}

		rigth = stack[--size];
		left = stack[--size];
	}
}

double TracingConcave::CalcMinDistanceToFacet(const Polygon &facet,
											  const Point3f &beamDir)
{
	double dist = FLT_MAX;
	const Point3f *pol = facet.arr;

	for (int i = 0; i < facet.size; ++i)
	{
		// measure dist
		Point3f point;
		ProjectPointToFacet(pol[i], -beamDir, beamDir, point);
		double newDist = sqrt(Norm(point - pol[i]));

		if (newDist < dist) // choose minimum with previews
		{
			dist = newDist;
		}
	}

	return dist;
}

void TracingConcave::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
										 std::vector<Beam> &/*outBeams*/)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		/// TODO: realize
	}
}

void TracingConcave::CutBeamByFacet(Paths &beamPolygon, int facetId,
									const Point3f &direction,
									const Point3f &facetNormal, Paths &result)
{
	Axis axis = m_clipper.GetSwapAxis(facetNormal);
	m_clipper.SwapCoords(axis, Axis::aZ, beamPolygon);

	Paths clip(1);
	{
		ProjectFacetToFacet(m_facets[facetId].polygon, direction, facetNormal,
							clip[0]); // проецируем грань на начальный пучок
		m_clipper.SwapCoords(axis, Axis::aZ, clip);
	}

	m_clipper.Difference(beamPolygon, clip, result);

	if (!result.empty())
	{
		m_clipper.HandleResultPaths(axis, result);
	}
}

double TracingConcave::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

//	Point3f normal;
//	Point3f p1 = beam.polygon[1] - beam.polygon[0];
//	Point3f p2 = beam.polygon[2] - beam.polygon[0];
//	CrossProduct(p1, p2, normal);

	// normal of last facet of beam
	Point3f normal = m_particle->facets[beam.facetId].ex_normal;

	double cosND = DotProduct(normal, beam.direction);
//	Normalize(normal);

//	if (fabs(normal.cz) < 0.5)
//	{
//		int minI = 0;

//		for (int i = 1; i < beam.shapeSize; ++i)
//		{
//			if (beam.shape[i].cz < beam.shape[minI].cz)
//			{
//				minI = i;
//			}
//		}
//	}
//	else if ...

//	if (dir < -FLT_EPSILON) /// нормаль посчитана неверно, разворачиваем её
//	{
//		normal.cx = -normal.cx;
//		normal.cy = -normal.cy;
//		normal.cz = -normal.cz;
//	}

	double e = fabs(cosND);

	if (e < eps)
	{
		return 0;
	}

	double square = m_clipper.AreaOfConcavePolygon(beam.polygon, normal);
	double n = Length(normal);
	return (e*square) / n;
}
