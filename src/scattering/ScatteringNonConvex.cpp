#include "ScatteringNonConvex.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>
#include <iostream>

#include "BigIntegerLibrary.hh"

//#ifdef _DEBUG // DEB
#include "Tracer.h"
//#endif

#define EPS_ORTO_FACET 0.0001

#ifdef _DEBUG
using namespace std;
ofstream trackMapFile("tracks_deb.dat", ios::out);
#endif

using namespace std;

ScatteringNonConvex::ScatteringNonConvex(Particle *particle, Light *incidentLight,
										 bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringNonConvex::ScatterLight(double beta, double gamma,
									   std::vector<Beam> &scaterredBeams)
{
	m_particle->Rotate(beta, gamma, 0);
	SplitLightToBeams();
	SplitBeams(scaterredBeams);
}

void ScatteringNonConvex::PushBeamsToTree(int facetId, const PolygonArray &polygons,
										  Beam &inBeam, Beam &outBeam)
{
	auto id = RecomputeTrackId(0, facetId);
	inBeam.SetTracingParams(facetId, 0, Location::In);
	outBeam.SetTracingParams(facetId, 0, Location::Out);
	inBeam.id = id;
	outBeam.id = id;

	for (unsigned j = 0; j < polygons.size; ++j)
	{
		Beam in = inBeam;
		Beam out = outBeam;
		const Polygon &pol = polygons.arr[j];
		in = pol;
		out = pol;
		Point3f p = pol.Center();
		in.opticalPath = 0;
		out.opticalPath = 0;
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
		in.AddOpticalPath(path);
		out.AddOpticalPath(path);
		m_beamTree[m_treeSize++] = in;
		m_beamTree[m_treeSize++] = out;
#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facetId, outBeam);
#endif
	}
}

void ScatteringNonConvex::SplitByFacet(const IntArray &facetIDs, int facetIndex)
{
	PolygonArray resPolygons;
	IntersectWithFacet(facetIDs, facetIndex, resPolygons);

	if (resPolygons.size != 0)
	{
		int id = facetIDs.arr[facetIndex];
		Beam inBeam, outBeam;
		SetIncidentBeamOpticalParams(id, inBeam, outBeam);
		PushBeamsToTree(id, resPolygons, inBeam, outBeam);
	}
}

void ScatteringNonConvex::SplitLightToBeams()
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
#endif
	m_treeSize = 0;

	IntArray facetIDs;
	SelectVisibleFacetsForLight(facetIDs);

	for (int i = 0; i < facetIDs.size; ++i)
	{
		SplitByFacet(facetIDs, i);
	}
}

void ScatteringNonConvex::SelectVisibleFacetsForLight(IntArray &facetIDs)
{
	FindVisibleFacetsForLight(facetIDs);
	SortFacets_faster(m_incidentDir, facetIDs);
}

void ScatteringNonConvex::IntersectWithFacet(const IntArray &facetIds, int prevFacetNum,
											 PolygonArray &resFacets)
{
	int id = facetIds.arr[prevFacetNum];

	if (prevFacetNum == 0 || m_facets[id].isVisibleOut)
	{
		resFacets.Push(m_facets[id]);
	}
	else // facet is probably shadowed by others
	{
		const Facet &facet = m_facets[id];
		const Point3f &normal = facet.ex_normal;

		CutPolygonByFacets(facet, facetIds, prevFacetNum, normal, normal,
						   m_incidentDir, resFacets);
	}
}

void ScatteringNonConvex::SelectVisibleFacets(const Beam &beam, IntArray &facetIDs)
{
	FindVisibleFacets(beam, facetIDs);

	Point3f dir = beam.direction;
	dir.d_param = m_facets[beam.lastFacetId].in_normal.d_param;
	SortFacets_faster(dir, facetIDs);
}

void ScatteringNonConvex::CutPolygonByFacets(const Polygon &pol,
											 const IntArray &facetIds, size_t size,
											 const Vector3f &polNormal,
											 const Vector3f &clipNormal,
											 const Vector3f &dir,
											 PolygonArray &pols)
{
	pols.Push(pol);

	// cut facet projections out of polygon one by one
	for (unsigned i = 0; i < size; ++i)
	{
		int id = facetIds.arr[i];
		m_polygonBuffer.Clear();

		while (pols.size != 0)
		{
			const Polygon &subj = pols.Pop();
			const Polygon &clip = m_facets[id];

			/// REF: объединить 2 первых аргумента и 2 вторых
			Difference(subj, polNormal, clip, clipNormal, dir, m_polygonBuffer);
		}

		if (m_polygonBuffer.size != 0)
		{
			for (unsigned i = 0; i < m_polygonBuffer.size; ++i)
			{
				pols.Push(m_polygonBuffer.arr[i]);
			}
		}
		else // beam has layed on the facet totally
		{
			break;
		}
	}
}

void ScatteringNonConvex::CutExternalBeam(const Beam &beam,
										  std::vector<Beam> &scaterredBeams)
{
	const Point3f &n1 = m_facets[beam.lastFacetId].ex_normal;
	const Point3f &n2 = m_facets[beam.lastFacetId].in_normal;

	IntArray facetIds;
	SelectVisibleFacets(beam, facetIds);

	PolygonArray resultBeams;
	CutPolygonByFacets(beam, facetIds, facetIds.size, n1, n2,
					   -beam.direction, resultBeams);

	Beam tmp = beam;
	tmp.opticalPath += m_splitting.ComputeOutgoingOpticalPath(tmp); // добираем оптический путь

	for (unsigned i = 0; i < resultBeams.size; ++i)
	{
		tmp.SetPolygon(resultBeams.arr[i]);
		scaterredBeams.push_back(tmp);
	}
}

void ScatteringNonConvex::SortFacets_faster(const Point3f &beamDir, IntArray &facetIDs)
{
	if (facetIDs.size == 0)
	{
		return;
	}

	int vertices[MAX_VERTEX_NUM];

	for (int i = 0; i < facetIDs.size; ++i)
	{
		const int &id = facetIDs.arr[i];
		vertices[i] = FindClosestVertex(m_facets[id], beamDir);
	}

	int left = 0;
	int rigth = facetIDs.size - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

	while (true)
	{
		int id = (left + rigth)/2;
		int i = left;
		int j = rigth;

		Point3f base = m_facets[facetIDs.arr[id]].arr[vertices[id]];

		while (i <= j)
		{
			Point3f vecB;
			double cosVN;

			do
			{
				vecB = base - m_facets[facetIDs.arr[i]].arr[vertices[i]];
				cosVN = DotProduct(vecB, beamDir);
				++i;
			}
			while (cosVN > FLT_EPSILON);
			--i;

			do
			{
				vecB = base - m_facets[facetIDs.arr[j]].arr[vertices[j]];
				cosVN = DotProduct(vecB, beamDir);
				--j;
			}
			while (cosVN < EPS_M_COS_90);
			++j;

			if (i <= j)	// exchange elems
			{
				float temp_d = vertices[i];
				vertices[i] = vertices[j];
				vertices[j] = temp_d;

				int temp_v = facetIDs.arr[i];
				facetIDs.arr[i] = facetIDs.arr[j];
				facetIDs.arr[j] = temp_v;

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

int ScatteringNonConvex::FindClosestVertex(const Polygon &facet, const Point3f &beamDir)
{
	int closest = 0;

	for (unsigned i = 1; i < facet.size; ++i)
	{
		Point3f v = facet.arr[closest] - facet.arr[i];
		double cosVD = DotProduct(v, beamDir);

		if (cosVD > FLT_EPSILON)
		{
			closest = i;
		}
	}

	return closest;
}

void ScatteringNonConvex::CutBeamByFacet(const Facet &facet, Beam &beam,
										 PolygonArray &result)
{
	const Location &loc = beam.location;
	const Facet &beamFacet = m_facets[beam.lastFacetId];

	if (loc == Location::In && beamFacet.isVisibleIn)
	{
		return;
	}

	const Point3f &facetNormal = (loc == Location::Out) ? -beamFacet.normal[loc]
														:  beamFacet.normal[loc];
	Difference(beam, beamFacet.normal[loc],
			   facet, facetNormal, -beam.direction, result);

	if (result.size == 0) // beam is totaly swallowed by facet
	{
		beam.size = 0;
	}
}

bool ScatteringNonConvex::IsOutgoingBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Out
			&& incidentBeam.size != 0); // OPT: replace each other
}

int ScatteringNonConvex::FindFacetId(int facetId, const IntArray &arr)
{
	int i = 0;

	while ((facetId == arr.arr[i]) && (i < arr.size))
	{
		++i;
	}

	if (i == arr.size)
	{
		i = -1;
	}

	return i;
}

void ScatteringNonConvex::SplitBeams(std::vector<Beam> &scaterredBeams)
{
#ifdef _DEBUG // DEB
	int count = 0;
#endif
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];

		if (!IsTerminalAct(beam)) // REF, OPT: перенести проверку во все места, где пучок закидывается в дерево, чтобы пучки заранее не закидывались в него
		{
			IntArray facetIds;
			SelectVisibleFacets(beam, facetIds);

			bool isDivided = false;

			for (unsigned i = 0; (i < facetIds.size) && !isDivided; ++i)// OPT: move this loop to SplitBeamByFacet
			{
				int facetId = facetIds.arr[i];
				Polygon intersection;
				Intersect(facetId, beam, intersection);

				if (intersection.size >= MIN_VERTEX_NUM)
				{
					isDivided = SplitBeamByFacet(intersection, facetId, beam);
				}
			}

			if (IsOutgoingBeam(beam))
			{	// посылаем обрезанный всеми гранями внешний пучок на сферу
				beam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
				scaterredBeams.push_back(beam);
			}
		}
		else if (beam.location == Location::Out)
		{
			CutExternalBeam(beam, scaterredBeams);
		}
	}
}

bool ScatteringNonConvex::SetOpticalBeamParams(const Facet &facet, const Beam &incidentBeam,
											   Beam &inBeam, Beam &outBeam)
{
	const Point3f &dir = incidentBeam.direction;
	const Point3f &normal = facet.ex_normal;

	bool hasOutBeam = true;
	m_splitting.ComputeCosA(dir, normal);

	if (m_splitting.IsNormalIncidence()) // normal incidence
	{
		m_splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);
	}
	else // regular incidence
	{
		Beam incBeam = incidentBeam;

		if (incidentBeam.location == Location::In)
		{
			m_splitting.ComputeSplittingParams(incidentBeam.direction, normal);
			ComputePolarisationParams(incBeam.direction, normal, incBeam);
//			incBeam.direction = -incBeam.direction;

			hasOutBeam = !m_splitting.IsCompleteReflection();

			if (hasOutBeam)
			{
				m_splitting.ComputeRegularBeamsParams(normal, incidentBeam,
													  inBeam, outBeam);
			}
			else // complete internal reflection incidence
			{
				m_splitting.ComputeCRBeamParams(normal, incidentBeam, inBeam);
			}
		}
		else // beam is external
		{
			inBeam.J = incidentBeam.J;
			m_splitting.ComputeCosA(dir, -normal);

			const Point3f &facetNormal = facet.in_normal;
			ComputePolarisationParams(-incBeam.direction, facetNormal, incBeam);
			m_splitting.ComputeRegularBeamParamsExternal(facetNormal, incBeam,
														 inBeam, outBeam);
//			if (m_isOpticalPath)
			{
				double path = m_splitting.ComputeSegmentOpticalPath(incidentBeam,
																	inBeam.Center());
				path += incidentBeam.opticalPath;
				inBeam.AddOpticalPath(path);
				outBeam.AddOpticalPath(path);
			}
		}
	}

	return hasOutBeam;
}

void ScatteringNonConvex::FindVisibleFacetsForLight(IntArray &facetIDs)
{
	for (int i = 0; i < m_particle->nFacets; ++i)
	{
		double cosA = DotProduct(m_incidentDir, m_facets[i].in_normal);

		if (cosA >= FLT_EPSILON) // beam incidents to this facet
		{
			facetIDs.Add(i);
		}
	}
}

bool ScatteringNonConvex::IsVisibleFacet(int facetID, const Beam &beam)
{
//	int loc = !beam.location;
	const Point3f &beamNormal = -m_facets[beam.lastFacetId].normal[!beam.location];

	const Point3f &facetCenter = m_facets[facetID].center;
	const Point3f &beamCenter = m_facets[beam.lastFacetId].center;
	Point3f vectorFromBeamToFacet = facetCenter - beamCenter;

	double cosBF = DotProduct(beamNormal, vectorFromBeamToFacet);
	return (cosBF >= EPS_ORTO_FACET);
}

void ScatteringNonConvex::FindVisibleFacets(const Beam &beam, IntArray &facetIds)
{
	int begin = 0;
	int end = m_particle->nFacets;

	if (m_particle->isAggregated && beam.location == Location::In)
	{
		m_particle->GetParticalFacetIdRangeByFacetId(beam.lastFacetId, begin, end);
	}

	for (int i = begin; i < end; ++i)
	{
		const Point3f &facetNormal = m_facets[i].normal[!beam.location];
		double cosFB = DotProduct(beam.direction, facetNormal);

		if (cosFB >= FLT_EPSILON) // beam incidents to this facet
		{
			if (IsVisibleFacet(i, beam))
			{	// facet is in front of begin of beam
				facetIds.Add(i);
			}
		}
	}
}

/// OPT: поменять все int и пр. параметры функций на ссылочные

/** TODO: придумать более надёжную сортировку по близости
 * (как вариант определять, что одна грань затеняют другую по мин. и макс.
 * удалённым вершинам, типа: "//" )
*/
void ScatteringNonConvex::SortFacets(const Point3f &beamDir, IntArray &facetIds)
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

double ScatteringNonConvex::CalcMinDistanceToFacet(const Polygon &facet,
											  const Point3f &beamDir)
{
	double dist = FLT_MAX;
	const Point3f *pol = facet.arr;
	Point3f point;
	Point3f dir = -beamDir;
	double dp = DotProduct(dir, beamDir);

	for (int i = 0; i < facet.size; ++i)
	{
		/// REF: заменить на сущ. фуyкцию ProjectPointToPlane
		// measure dist
		double t = DotProduct(pol[i], beamDir);
		t = t + beamDir.d_param;
		t = t/dp;
		point = pol[i] - (dir * t);
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

void ScatteringNonConvex::PushBeamPartsToTree(const Beam &beam,
											  const PolygonArray &parts)
{
	Beam tmp = beam; // OPT: try to replace 'tmp' to 'beam'

	for (unsigned i = 0; i < parts.size; ++i)
	{
		tmp = parts.arr[i];
		assert(m_treeSize < MAX_BEAM_REFL_NUM);
		m_beamTree[m_treeSize++] = tmp;
	}
}

void ScatteringNonConvex::PushBeamToTree(Beam &beam, const Beam &oldBeam,
										 const BigInteger &newId, int facetId,
										 Location loc)
{
	beam.id = newId;
	beam.locations = oldBeam.locations;
	Scattering::PushBeamToTree(beam, facetId, oldBeam.nActs+1, loc);
}

bool ScatteringNonConvex::SplitBeamByFacet(const Polygon &intersection,
										   int facetId, Beam &beam)
{
#ifdef _DEBUG // DEB
if (beam.lastFacetId==0 && facetId==6)
//if (beam.trackId.toLong()==9633 && facetID==0)
	int f =0;
#endif
	Facet &facet = m_facets[facetId];

	Beam inBeam, outBeam;
	inBeam.SetPolygon(intersection);
	outBeam.SetPolygon(intersection);

	bool hasOutBeam = SetOpticalBeamParams(facet, beam, inBeam, outBeam);

	auto newId = RecomputeTrackId(beam.id, facetId);

	if (hasOutBeam)
	{
		PushBeamToTree(outBeam, beam, newId, facetId, Location::Out);
	}

	PushBeamToTree(inBeam, beam, newId, facetId, Location::In);

	m_polygonBuffer.Clear();
	CutBeamByFacet(facet, beam, m_polygonBuffer);

	bool isDivided = m_polygonBuffer.size > CLIP_RESULT_SINGLE;

	if (isDivided)
	{	// beam had divided by facet
		PushBeamPartsToTree(beam, m_polygonBuffer);
		beam.size = 0;
	}
	else if (m_polygonBuffer.size == CLIP_RESULT_SINGLE)
	{
		beam = m_polygonBuffer.arr[0];
	}

	return isDivided;
}

void ScatteringNonConvex::PushBeamsToBuffer(int facetID, const Beam &beam, bool hasOutBeam,
									   Beam &inBeam, Beam &outBeam,
									   std::vector<Beam> &passed)
{
	inBeam.id = beam.id;

	if (hasOutBeam)
	{
		outBeam.id = beam.id;
		outBeam.SetTracingParams(facetID, beam.nActs+1, Location::Out);
		outBeam.id = RecomputeTrackId(outBeam.id, outBeam.lastFacetId);
		passed.push_back(outBeam);
	}

	inBeam.SetTracingParams(facetID, beam.nActs+1, Location::In);
	inBeam.id = RecomputeTrackId(outBeam.id, outBeam.lastFacetId);
	passed.push_back(inBeam);
}

void ScatteringNonConvex::ScatterLight(double beta, double gamma,
										 const std::vector<std::vector<int>> &tracks,
										 std::vector<Beam> &scaterredBeams)
{
//	m_particle->Rotate(beta, gamma, 0);

//	for (const std::vector<int> &track : tracks)
//	{
//		int facetID = track.at(0);

//		bool isIncident;
//		TraceFirstBeamFixedFacet(facetID, isIncident);

//		if (!isIncident)
//		{
//			continue;
//		}

//		for (size_t i = 1; i < track.size(); ++i)
//		{
//			int facetID = track.at(i);

//			std::vector<Beam> buffer; // для прошедших пучков (не дублированных)

//			while (m_treeSize != 0)
//			{
//				Beam beam = m_beamTree[--m_treeSize];

//				IntArray facetIDs;
//				SelectVisibleFacets(beam, facetIDs);
//				int index = FindFacetId(facetID, facetIDs);

//				if (index != -1)
//				{
//					bool isDivided;
//					SplitBeamByFacet(beam, facetID, isDivided);

//					Polygon intersected;
//					Intersect(facetID, beam, intersected);

//					if (intersected.size >= MIN_VERTEX_NUM)
//					{
//						Beam inBeam, outBeam;
//						inBeam.SetPolygon(intersected);
//						outBeam.SetPolygon(intersected);

//						Facet &facet = m_facets[facetID];
//						bool hasOutBeam = SetOpticalBeamParams(facet, beam,
//															   inBeam, outBeam);
//						PushBeamsToBuffer(facetID, beam, hasOutBeam,
//										  inBeam, outBeam, buffer);
//					}
//				}
//			}

//			if (buffer.empty())
//			{
//				break;
//			}

//			for (const Beam &b : buffer)
//			{	// добавляем прошедшие пучки
//				assert(m_treeSize < MAX_BEAM_REFL_NUM);
//				m_beamTree[m_treeSize++] = b;
//			}
//		}

//		while (m_treeSize != 0)
//		{
//			scaterredBeams.push_back(m_beamTree[--m_treeSize]);
//		}
//	}
}

void ScatteringNonConvex::TraceFirstBeamFixedFacet(int facetID, bool &isIncident)
{
	IntArray facetIDs;
	isIncident = false;

	SelectVisibleFacetsForLight(facetIDs);
	int index = FindFacetId(facetID, facetIDs);

	if (index != -1)
	{
		SplitByFacet(facetIDs, index);
		isIncident = true;
	}
}
