#include "TracingConcave.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>

#ifdef _DEBUG // DEB
#include <iostream>
#endif

#define EPS_ORTO_FACET 0.0001 //TODO: подобрать норм значение

#ifdef _TRACK_ALLOW
std::ofstream trackMapFile("tracks.dat", std::ios::out);
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
	TraceSecondaryBeams(scaterredBeams);

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

	m_treeSize = 0;
	m_lightSurfaceArea = 0;

	IntArray facetIDs;
	SelectVisibleFacetsForWavefront(facetIDs);

	for (int i = 0; i < facetIDs.size; ++i)
	{
		PolygonArray resFacets;
		bool hasIntersection;
		IntersectWithFacet(facetIDs, i, resFacets, hasIntersection);

		if (hasIntersection)
		{
			Beam inBeam, outBeam;

			// set optical params of beam
			int facetId = facetIDs.arr[i];
			SetFirstBeamOpticalParams(facetId, inBeam, outBeam);

#ifdef _DEBUG // DEB
			double eee = 0;
#endif
			for (int j = 0; j < resFacets.size; ++j)
			{
#ifdef _DEBUG // DEB
				eee += AreaOfPolygon(resFacets[j]);
#endif
				// set geometry of beam
				 inBeam.SetPolygon(resFacets.arr[j]);
				outBeam.SetPolygon(resFacets.arr[j]);

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

void TracingConcave::SelectVisibleFacetsForWavefront(IntArray &facetIds)
{
	FindVisibleFacetsForWavefront(facetIds);
	SortFacets(m_waveFront.direction, facetIds);
}

void TracingConcave::IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
										PolygonArray &resFacets, bool &hasIntersection)
{
	resFacets.size = 0;
	hasIntersection = true;
	int facetId = facetIds.arr[prevFacetNum];

	if (prevFacetNum == 0 || m_particle->IsUnshadowedExternal(facetId))
	{
		Polygon beamPolygon;
		SetPolygonByFacet(facetId, beamPolygon);
		resFacets.arr[resFacets.size++] = beamPolygon;
	}
	else // facet is probably shadowed by others
	{
		CutShadowsFromFacet(facetId, facetIds, prevFacetNum, m_waveFront,
							resFacets);

		if (resFacets.size == 0)
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
	const Point3f &normal = m_facets[beam.facetId].in_normal;
//	const Point3f &facetNormal = m_facets[beam.facetId].ex_normal;
//	const Point3f &facetNormal = (beam.location == Location::Outside) ?  normal
//																	  : -normal;
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
			const Polygon &clip = m_facets[id].polygon;
			const Polygon &subj = resultBeams[--resSize];
//			Difference(subj, normal, clip, facetNormal, -beam.direction,
//					   diffFacets, diffSize);
			Difference(subj, normal, clip, normal, -beam.direction,
					   diffFacets, diffSize);
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

void TracingConcave::CutBeamByFacet(int facetId, Beam &beam, bool &isDivided)
{
	isDivided = false;
	const Location &loc = beam.location;

	if (loc == Location::Inside && !m_particle->IsShadowedInternal(beam.facetId))
	{
		return;
	}

	const Facet &beamFacet = m_facets[beam.facetId];
	const Point3f &beamNormal = beamFacet.normal[loc];
	const Point3f &facetNormal = (loc == Location::Outside) ? -beamNormal
															:  beamNormal;
	Polygon resultBeams[MAX_VERTEX_NUM];
	int resultSize = 0;
	Difference(beam.polygon, beamNormal, m_facets[facetId].polygon, facetNormal,
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
	outBeam.track = incidentBeam.track;
	inBeam.track = incidentBeam.track;
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
	const Point3f &beamNormal = -m_facets[beam.facetId].normal[!beam.location];

	const Point3f &facetCenter = m_particle->centers[facetID];
	const Point3f &beamCenter = m_particle->centers[beam.facetId];
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

void TracingConcave::CutShadowsFromFacet(int facetId, const IntArray &facetIds,
										 int prevFacetNum, const Beam &beam,
										 PolygonArray &resFacets)
{
	const Facet &facet = m_facets[facetId];
	const Point3f &facetNormal = facet.normal[Location::Outside];
	const Point3f &clipNormal = facetNormal;

	resFacets.arr[resFacets.size++] = facet.polygon;

	for (int i = 0; i < prevFacetNum; ++i)
	{
		int id = facetIds.arr[i];

		Polygon diffFacets[MAX_POLYGON_NUM];
		int diffSize = 0;

		while (resFacets.size != 0)
		{
			const Polygon &clip = m_facets[id].polygon;
			const Polygon &subj = resFacets.arr[--resFacets.size];
//			const Point3f &clipNormal = m_facets[id].normal[(int)beam.location];

			Difference(subj, facetNormal, clip, clipNormal, -beam.direction,
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

	Point3f normal = m_facets[beam.facetId].ex_normal; // normal of last facet of beam
	double cosFB = DotProduct(normal, beam.direction);
	double e = fabs(cosFB);

	if (e < eps)
	{
		return 0;
	}

	double square = AreaOfPolygon(beam.polygon);
	double n = Length(normal);
	return (e*square) / n;
}
