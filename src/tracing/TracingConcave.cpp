#include "TracingConcave.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>

#ifdef _DEBUG // DEB
#include <iostream>
#endif

#define EPS_ORTO_FACET 0.0001

#ifdef _TRACK_ALLOW
std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
#endif

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	Point3f point = m_waveFront.direction * m_particle->GetHalfHeight()*2;
	m_waveFront.direction.d_param = DotProduct(point, m_waveFront.direction);
	m_waveFront.location = Location::Outside;

	m_isArea = false;
}

void TracingConcave::SplitBeamByParticle(std::vector<Beam> &scaterredBeams)
{
	TraceFirstBeam();

#ifdef _DEBUG // DEB
	double rrr = 0;
	for (int i = 0; i < m_treeSize; ++i)
	{
//		if (b.track.size() == 1 &&
//				b.track[0] == 1)
			rrr += m_beamTree[i].polygon.Area();
	}
	int fgfg = 0;
#endif

	TraceSecondaryBeams(scaterredBeams);

#ifdef _DEBUG // DEB
	double fff = 0;
	for (const Beam &b : scaterredBeams)
	{
//		if (b.track.size() == 1 &&
//				b.track[0] == 1)
			fff += b.polygon.Area();
	}
	int ggg = 0;
#endif
}

void TracingConcave::PushBeamsToTree(int facetID, const PolygonArray &polygons,
									 Beam &inBeam, Beam &outBeam)
{
	for (int j = 0; j < polygons.size; ++j)
	{
		// set geometry of beam
		 inBeam.SetPolygon(polygons.arr[j]);
		outBeam.SetPolygon(polygons.arr[j]);

		PushBeamToTree( inBeam, facetID, 0, Location::Inside );
		PushBeamToTree(outBeam, facetID, 0, Location::Outside);

		if (m_isArea)
		{
			CalcLigthSurfaceArea(facetID, outBeam);
		}
	}
}

void TracingConcave::TraceFirstBeam()
{
#ifdef _TRACK_OUTPUT
	trackMapFile << "0 lvl: ";
#endif

	m_treeSize = 0;
	m_lightSurfaceArea = 0;

	IntArray facetIDs;
	SelectVisibleFacetsForWavefront(facetIDs);

	for (int i = 0; i < facetIDs.size; ++i)
	{
		PolygonArray resPolygons;
		IntersectWithFacet(facetIDs, i, resPolygons);

		if (resPolygons.size != 0)
		{
			int id = facetIDs.arr[i];
			Beam inBeam, outBeam;
			SetBeamOpticalParams(id, inBeam, outBeam);
			PushBeamsToTree(id, resPolygons, inBeam, outBeam);
		}
	}
}

void TracingConcave::SelectVisibleFacetsForWavefront(IntArray &facetIds)
{
	FindVisibleFacetsForWavefront(facetIds);
	SortFacets(m_waveFront.direction, facetIds);
}

void TracingConcave::IntersectWithFacet(const IntArray &facetIDs, int prevFacetNum,
										PolygonArray &resFacets)
{
	int id = facetIDs.arr[prevFacetNum];

	if (prevFacetNum == 0 || m_particle->IsUnshadowedExternal(id))
	{
		resFacets.arr[resFacets.size++] = m_facets[id].polygon;
	}
	else // facet is probably shadowed by others
	{
		CutFacetByShadows(id, facetIDs, prevFacetNum, resFacets);
	}
}

void TracingConcave::SelectVisibleFacets(const Beam &beam, IntArray &facetIds)
{
	FindVisibleFacets(beam, facetIds);

	Point3f dir = beam.direction;
	dir.d_param = m_facets[beam.facetID].in_normal.d_param;
	SortFacets(dir, facetIds);
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
{
	const Point3f &normal = m_facets[beam.facetID].ex_normal;
	const Point3f &normal1 = m_facets[beam.facetID].in_normal;

	IntArray facetIds;
	SelectVisibleFacets(beam, facetIds);

	Polygon resultBeams[MAX_VERTEX_NUM];
	int resSize = 0;
	resultBeams[resSize++] = beam.polygon;

	// cut facet projections out of beam one by one
	for (int i = 0; i < facetIds.size; ++i)
	{
		int id = facetIds.arr[i];

		Polygon diffFacets[MAX_POLYGON_NUM];
		int diffSize = 0;

		while (resSize != 0)
		{
			Difference(resultBeams[--resSize], normal, m_facets[id].polygon,
					normal1, -beam.direction, diffFacets, diffSize);
		}

		if (diffSize != 0)
		{
			for (int i = 0; i < diffSize; ++i)
			{
				resultBeams[resSize++] = diffFacets[i];
			}
		}
		else // beam is totaly swallowed by facet
		{
			break;
		}
	}

	Beam tmp = beam;

	for (int i = 0; i < resSize; ++i)
	{
		tmp.SetPolygon(resultBeams[i]);
		scatteredBeams.push_back(tmp);
	}
}

#ifdef _TRACK_ALLOW
//void TracingConcave::AddToTrack(Beam &beam, int facetId)
//{
//	std::vector<int> &tr = beam.track;
//	int size = tr.size();

//	if (size == 0 || (size > 0 && facetId != tr.at(size-1)))
//		tr.push_back(facetId);
//}

//void TracingConcave::PrintTrack(const Beam &beam, int facetId)
//{
//	trackMapFile << "(" << beam.track.at(0);

//	for (int i = 1; i < beam.track.size(); ++i)
//	{
//		trackMapFile << ", ";
//		trackMapFile << beam.track.at(i);
//	}

//	if (facetId != beam.track.back())
//	{
//		trackMapFile << ", " << facetId;
//	}

//	trackMapFile << "), ";
//	trackMapFile.flush();
//}
#endif

void TracingConcave::CutBeamByFacet(int facetId, Beam &beam, bool &isDivided)
{
	isDivided = false;
	const Location &loc = beam.location;

	if (loc == Location::Inside && !m_particle->IsShadowedInternal(beam.facetID))
	{
		return;
	}

	const Facet &beamFacet = m_facets[beam.facetID];
	const Point3f &facetNormal = (loc == Location::Outside) ? -beamFacet.normal[loc]
															:  beamFacet.normal[loc];
	Polygon resultBeams[MAX_VERTEX_NUM];
	int resultSize = 0;
	Difference(beam.polygon, beamFacet.normal[loc],
			   m_facets[facetId].polygon, facetNormal, -beam.direction,
			   resultBeams, resultSize);

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
			PushBeamToTree(tmp);
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

void TracingConcave::PushBeamsToTree(int level, int facetID, bool hasOutBeam,
									 Beam &inBeam, Beam &outBeam)
{
	if (hasOutBeam)
	{
		PushBeamToTree(outBeam, facetID, level, Location::Outside);
	}

#ifdef _TRACK_ALLOW
	inBeam.id = outBeam.id;
#ifdef _TRACK_OUTPUT
	trackMapFile << "[in] ";
#endif
#endif
	PushBeamToTree(inBeam, facetID, level, Location::Inside);
}

void TracingConcave::TraceSecondaryBeams(std::vector<Beam> &scaterredBeams)
{
	while (m_treeSize != 0)
	{
		Beam incidentBeam = m_beamTree[--m_treeSize];
		const Location loc = incidentBeam.location;

		if (IsTerminalBeam(incidentBeam))
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
				CutBeamByFacet(facetId, incidentBeam, isDivided);

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
	outBeam.id = incidentBeam.id;
	inBeam.id = incidentBeam.id;
#endif
	const Point3f &incidentDir = incidentBeam.direction;
	const Point3f &normal = m_facets[facetId].ex_normal;

	hasOutBeam = true;

	double cosIN = DotProduct(incidentDir, normal);

	if (cosIN >= EPS_COS_00) // normal incidence
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
			inBeam.J = incidentBeam.J;
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

void TracingConcave::FindVisibleFacetsForWavefront(IntArray &facetIds)
{
	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIN = DotProduct(m_waveFront.direction, m_facets[i].in_normal);

		if (cosIN >= EPS_COS_90) // beam incidents to this facet
		{
			facetIds.arr[facetIds.size++] = i;
		}
	}
}

bool TracingConcave::IsVisibleFacet(int facetID, const Beam &beam)
{
//	int loc = !beam.location;
	const Point3f &beamNormal = -m_facets[beam.facetID].normal[!beam.location];

	const Point3f &facetCenter = m_particle->centers[facetID];
	const Point3f &beamCenter = m_particle->centers[beam.facetID];
	Point3f vectorFromBeamToFacet = facetCenter - beamCenter;

	double cosBF = DotProduct(beamNormal, vectorFromBeamToFacet);
	return (cosBF >= EPS_ORTO_FACET);
}

void TracingConcave::FindVisibleFacets(const Beam &beam, IntArray &facetIds)
{
	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		const Point3f &facetNormal = m_facets[i].normal[!beam.location];
		double cosFB = DotProduct(beam.direction, facetNormal);

		if (cosFB >= EPS_COS_90) // beam incidents to this facet
		{
			if (IsVisibleFacet(i, beam))
			{	// facet is in front of begin of beam
				facetIds.arr[facetIds.size++] = i;
			}
		}
	}
}

void TracingConcave::CutFacetByShadows(int facetID, const IntArray &shadowFacetIDs,
									   int prevFacetNum, PolygonArray &resFacets)
{
	const Facet &facet = m_facets[facetID];
	const Point3f &normal = facet.normal[Location::Outside];
	resFacets.arr[resFacets.size++] = facet.polygon;

	for (int i = 0; i < prevFacetNum; ++i)
	{
		int id = shadowFacetIDs.arr[i];
		Polygon diffFacets[MAX_POLYGON_NUM];
		int diffSize = 0;

		while (resFacets.size != 0)
		{
			const Polygon &clip = m_facets[id].polygon;
			const Polygon &subj = resFacets.arr[--resFacets.size];
			Difference(subj, normal, clip, normal, -m_waveFront.direction,
					   diffFacets, diffSize);
		}

		if (diffSize != 0)
		{
			for (int i = 0; i < diffSize; ++i)
			{
				resFacets.arr[resFacets.size++] = diffFacets[i];
			}
		}
		else // beam is totaly swallowed by facet
		{
			break;
		}
	}
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
										 std::vector<Beam> &/*scaterredBeams*/)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		/// TODO: realize
	}
}

double TracingConcave::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f normal = m_facets[beam.facetID].ex_normal; // normal of last facet of beam
	double cosFB = DotProduct(normal, beam.direction);
	double e = fabs(cosFB);

	if (e < eps)
	{
		return 0;
	}

	double square = beam.polygon.Area();
	double n = Length(normal);
	return (e*square) / n;
}
