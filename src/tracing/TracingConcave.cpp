#include "TracingConcave.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>

#ifdef _DEBUG // DEB
#include <iostream>
#endif

#define EPS_ORTO_FACET 0.0001

#ifdef _TRACK_ALLOW
//std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
#endif

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	Point3f point = m_incidentDir * m_particle->GetMainSize();
	m_incidentDir.d_param = DotProduct(point, m_incidentDir);
}

void TracingConcave::SplitBeamByParticle(double beta, double gamma,
										 std::vector<Beam> &scaterredBeams)
{
	m_particle->Rotate(beta, gamma, 0);
	TraceFirstBeam();
	TraceSecondaryBeams(scaterredBeams);
}

void TracingConcave::PushBeamsToTree(int facetID, const PolygonArray &polygons,
									 Beam &inBeam, Beam &outBeam)
{
	for (int j = 0; j < polygons.size; ++j)
	{
		long long inID = inBeam.id;
		long long outID = outBeam.id;

		// set geometry of beam
		 inBeam.SetPolygon(polygons.arr[j]);
		outBeam.SetPolygon(polygons.arr[j]);

		// REF: дублирует одиночный вызов в пред. ф-ции
		CalcOpticalPath_initial(inBeam, outBeam);

		PushBeamToTree( inBeam, facetID, 0, Location::In );
		PushBeamToTree(outBeam, facetID, 0, Location::Out);

#ifdef _CHECK_ENERGY_BALANCE
		CalcFacetEnergy(facetID, outBeam);
#endif
		inBeam.id = inID;
		outBeam.id = outID;
	}
}

void TracingConcave::TraceByFacet(const IntArray &facetIDs, int facetIndex)
{
#ifdef _DEBUG // DEB
	if (facetIndex == 25)
		int fff = 0;
#endif
	PolygonArray resPolygons;
	IntersectWithFacet(facetIDs, facetIndex, resPolygons);

	if (resPolygons.size != 0)
	{
		int id = facetIDs.arr[facetIndex];
		Beam inBeam, outBeam;
		SetBeamOpticalParams(id, inBeam, outBeam);
		PushBeamsToTree(id, resPolygons, inBeam, outBeam);
	}
}

void TracingConcave::TraceFirstBeam()
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incommingEnergy = 0;
#endif
	m_treeSize = 0;

	IntArray facetIDs;
	SelectVisibleFacetsForWavefront(facetIDs);

	for (int i = 0; i < facetIDs.size; ++i)
	{
		TraceByFacet(facetIDs, i);
	}
}

void TracingConcave::SelectVisibleFacetsForWavefront(IntArray &facetIDs)
{
	FindVisibleFacetsForWavefront(facetIDs);
	SortFacets(m_incidentDir, facetIDs);
}

void TracingConcave::IntersectWithFacet(const IntArray &facetIDs, int prevFacetNum,
										PolygonArray &resFacets)
{
	int id = facetIDs.arr[prevFacetNum];

	if (prevFacetNum == 0 || m_facets[id].isVisibleOut)
	{
		resFacets.arr[resFacets.size++] = m_facets[id];
	}
	else // facet is probably shadowed by others
	{
		CutFacetByShadows(id, facetIDs, prevFacetNum, resFacets);
	}
}

void TracingConcave::SelectVisibleFacets(const Beam &beam, IntArray &facetIDs)
{
	FindVisibleFacets(beam, facetIDs);

	Point3f dir = beam.direction;
	dir.d_param = m_facets[beam.lastFacetID].in_normal.d_param;
	SortFacets(dir, facetIDs);
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
{
	const Point3f &normal = m_facets[beam.lastFacetID].ex_normal;
	const Point3f &normal1 = m_facets[beam.lastFacetID].in_normal;

	IntArray facetIds;
	SelectVisibleFacets(beam, facetIds);

	Polygon resultBeams[MAX_VERTEX_NUM];
	int resSize = 0;
	resultBeams[resSize++] = beam;

	// cut facet projections out of beam one by one
	for (int i = 0; i < facetIds.size; ++i)
	{
		int id = facetIds.arr[i];

		Polygon diffFacets[MAX_POLYGON_NUM];
		int diffSize = 0;

		while (resSize != 0)
		{
			Difference(resultBeams[--resSize], normal, m_facets[id],
					normal1, -beam.direction, diffFacets, diffSize);
		}
#ifdef _DEBUG // DEB
		if (diffSize >= MAX_POLYGON_NUM)
			int fff = 0;
#endif
assert(diffSize < MAX_POLYGON_NUM); // DEB

		if (diffSize == 0) // beam is totaly swallowed by facet
		{
			break;
		}

		for (int i = 0; i < diffSize; ++i)
		{
			resultBeams[resSize++] = diffFacets[i];
assert((resSize < MAX_VERTEX_NUM) && resSize); // DEB
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
#endif

void TracingConcave::CutBeamByFacet(int facetID, Beam &beam, bool &isDivided)
{
	isDivided = false;
	const Location &loc = beam.location;
	const Facet &beamFacet = m_facets[beam.lastFacetID];

	if (loc == Location::In && beamFacet.isVisibleIn)
	{
		return;
	}

	const Point3f &facetNormal = (loc == Location::Out) ? -beamFacet.normal[loc]
														:  beamFacet.normal[loc];
	Polygon resultBeams[MAX_VERTEX_NUM];
	int resultSize = 0;
	Difference(beam, beamFacet.normal[loc],
			   m_facets[facetID], facetNormal, -beam.direction,
			   resultBeams, resultSize);

	if (resultSize == 0) // beam is totaly swallowed by facet
	{
		beam.size = 0;
	}
	else if (resultSize == CLIP_RESULT_SINGLE)
	{
		beam = resultBeams[0];
	}
	else // beam had divided by facet
	{
		Beam tmp = beam;

		for (int i = 0; i < resultSize; ++i)
		{
			tmp = resultBeams[i];
			m_beamTree[m_treeSize++] = tmp;
		}

		isDivided = true;
		beam.size = 0;
	}
}

bool TracingConcave::isExternalNonEmptyBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Out
			&& incidentBeam.size != 0);
}

int TracingConcave::FindFacetID(int facetID, const IntArray &arr)
{
	int i = 0;

	while ((facetID == arr.arr[i]) && (i < arr.size))
	{
		++i;
	}

	if (i == arr.size)
	{
		i = -1;
	}

	return i;
}

void TracingConcave::PushBeamsToTree(const Beam &beam, int facetID, bool hasOutBeam,
									 Beam &inBeam, Beam &outBeam)
{
assert((m_treeSize < MAX_BEAM_REFL_NUM)); // DEB
#ifdef _TRACK_ALLOW
	inBeam.id = beam.id;
#endif

	if (hasOutBeam)
	{
#ifdef _TRACK_ALLOW
		outBeam.id = beam.id;
#endif
		PushBeamToTree(outBeam, facetID, beam.level+1, Location::Out);
	}

	PushBeamToTree(inBeam, facetID, beam.level+1, Location::In);
}

void TracingConcave::TraceSecondaryBeams(std::vector<Beam> &scaterredBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];
#ifdef _DEBUG // DEB
		if (beam.lastFacetID == 13)
			int fff = 0;
#endif
		if (IsTerminalBeam(beam))
		{
			if (beam.location == Location::Out)
			{
				CatchExternalBeam(beam, scaterredBeams);
			}

			continue;
		}

		IntArray facetIds;
		SelectVisibleFacets(beam, facetIds);

		for (int i = 0; i < facetIds.size; ++i)
		{
			int facetID = facetIds.arr[i];

			bool isDivided;
			TraceSecondaryBeamByFacet(beam, facetID, isDivided);

			if (isDivided)
			{
				break;
			}
		}

		if (isExternalNonEmptyBeam(beam))
		{	// посылаем обрезанный всеми гранями внешний пучок на сферу
			scaterredBeams.push_back(beam);
		}
	}
}

void TracingConcave::SetOpticalBeamParams(int facetID, const Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam,
										  bool &hasOutBeam)
{
	const Point3f &incidentDir = incidentBeam.direction;
	const Point3f &normal = m_facets[facetID].ex_normal;

	hasOutBeam = true;
	double cosIN = DotProduct(incidentDir, normal);

	if (cosIN >= EPS_COS_00) // normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);
	}
	else // slopping incidence
	{
		if (incidentBeam.location == Location::In)
		{
			Beam incBeam = incidentBeam;
			SetSloppingIncidenceBeamParams(cosIN, normal, incBeam,
										   inBeam, outBeam, hasOutBeam);
		}
		else // beam is external
		{
			inBeam.J = incidentBeam.J;
			double cosI = DotProduct(-normal, incidentDir);

			SetSloppingBeamParams_initial(incidentDir, cosI, facetID,
										  inBeam, outBeam);
			if (m_isOpticalPath)
			{
				CalcOpticalPathInternal(cosIN, incidentBeam, inBeam, outBeam);
			}
		}
	}
}

void TracingConcave::FindVisibleFacetsForWavefront(IntArray &facetIDs)
{
	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIN = DotProduct(m_incidentDir, m_facets[i].in_normal);

		if (cosIN >= EPS_COS_90) // beam incidents to this facet
		{
			facetIDs.arr[facetIDs.size++] = i;
		}
	}
}

bool TracingConcave::IsVisibleFacet(int facetID, const Beam &beam)
{
//	int loc = !beam.location;
	const Point3f &beamNormal = -m_facets[beam.lastFacetID].normal[!beam.location];

	const Point3f &facetCenter = m_facets[facetID].center;
	const Point3f &beamCenter = m_facets[beam.lastFacetID].center;
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
	const Point3f &normal = facet.normal[Location::Out];
	resFacets.arr[resFacets.size++] = facet;

	for (int i = 0; i < prevFacetNum; ++i)
	{
		int id = shadowFacetIDs.arr[i];
		Polygon diffFacets[MAX_POLYGON_NUM]; // REF: заменить на структуру с size
		int diffSize = 0;
#ifdef _DEBUG // DEB
		if (i == 24)
			int fff = 0;
#endif
		while (resFacets.size != 0)
		{
			const Polygon &clip = m_facets[id];
			const Polygon &subj = resFacets.arr[--resFacets.size];
			Difference(subj, normal, clip, normal, /*m_facets[id].in_normal*/m_incidentDir,
					   diffFacets, diffSize);
		}

		if (diffSize == 0) // beam is totaly swallowed by facet
		{
			break;
		}

		for (int i = 0; i < diffSize; ++i)
		{
			resFacets.arr[resFacets.size++] = diffFacets[i];
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

// OPT: поменять все int и пр. параметры функций на ссылочные

/* TODO: придумать более надёжную сортировку по близости
 * (как вариант определять, что одна грань затеняют другую по мин. и макс.
 * удалённым вершинам, типа: "//" )*/
void TracingConcave::SortFacets(const Point3f &beamDir, IntArray &facetIds)
{
	float distances[MAX_VERTEX_NUM];

	for (int i = 0; i < facetIds.size; ++i)
	{
		const int &id = facetIds.arr[i];
		distances[i] = CalcMinDistanceToFacet(m_facets[id], beamDir);
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

/* TODO
 * Разобраться с параметром 'n' (кол-во вн. столкновений)
 * при заданных траекториях, возможно он не нужен т.к. заранее известен путь
 */

void TracingConcave::TraceSecondaryBeamByFacet(Beam &beam, int facetID,
											   bool &isDivided)
{
	isDivided = false;
	Polygon intersected;
	bool hasIntersection = Intersect(facetID, beam, intersected);

	if (hasIntersection)
	{
		Beam inBeam, outBeam;
		inBeam.SetPolygon(intersected);
		outBeam.SetPolygon(intersected);

		bool hasOutBeam;
		SetOpticalBeamParams(facetID, beam, inBeam, outBeam, hasOutBeam);
		PushBeamsToTree(beam, facetID, hasOutBeam, inBeam, outBeam);
		CutBeamByFacet(facetID, beam, isDivided);
	}
}

void TracingConcave::PushBeamsToBuffer(int facetID, const Beam &beam, bool hasOutBeam,
									   Beam &inBeam, Beam &outBeam,
									   std::vector<Beam> &passed)
{
	inBeam.id = beam.id;

	if (hasOutBeam)
	{
		outBeam.id = beam.id;
		outBeam.SetTracingParams(facetID, beam.level+1, Location::Out);
		SetBeamID(outBeam);
		passed.push_back(outBeam);
	}

	outBeam.SetTracingParams(facetID, beam.level+1, Location::In);
	SetBeamID(inBeam);
	passed.push_back(inBeam);
}

void TracingConcave::SplitBeamByParticle(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
										 std::vector<Beam> &scaterredBeams)
{
	m_particle->Rotate(beta, gamma, 0);

	for (const std::vector<int> &track : tracks)
	{
		int facetID = track.at(0);

		bool isIncident;
		TraceFirstBeamFixedFacet(facetID, isIncident);

		if (!isIncident)
		{
			continue;
		}

		for (int i = 0; i < track.size(); ++i)
		{
			int facetID = track.at(i);

			std::vector<Beam> buffer; // для прошедших пучков (не дублированных)

			while (m_treeSize != 0)
			{
				Beam beam = m_beamTree[--m_treeSize];

				IntArray facetIDs;
				SelectVisibleFacets(beam, facetIDs);
				int index = FindFacetID(facetID, facetIDs);

				if (index != -1)
				{
					bool isDivided;
					TraceSecondaryBeamByFacet(beam, facetID, isDivided);

					Polygon intersected;
					bool hasIntersection = Intersect(facetID, beam, intersected);

					if (hasIntersection)
					{
						Beam inBeam, outBeam;
						inBeam.SetPolygon(intersected);
						outBeam.SetPolygon(intersected);

						bool hasOutBeam;
						SetOpticalBeamParams(facetID, beam, inBeam, outBeam, hasOutBeam);
						PushBeamsToBuffer(facetID, beam, hasOutBeam, inBeam, outBeam, buffer);
						CutBeamByFacet(facetID, beam, isDivided);
					}
				}
			}

			if (buffer.empty())
			{
				break;
			}

			for (const Beam &b : buffer)
			{	// добавляем прошедшие пучки
				m_beamTree[m_treeSize++] = b;
			}
		}

		while (m_treeSize != 0)
		{
			scaterredBeams.push_back(m_beamTree[--m_treeSize]);
		}
	}
}

void TracingConcave::TraceFirstBeamFixedFacet(int facetID, bool &isIncident)
{
	IntArray facetIDs;
	isIncident = false;

	SelectVisibleFacetsForWavefront(facetIDs);
	int index = FindFacetID(facetID, facetIDs);

	if (index != -1)
	{
		TraceByFacet(facetIDs, index);
		isIncident = true;
	}
}
