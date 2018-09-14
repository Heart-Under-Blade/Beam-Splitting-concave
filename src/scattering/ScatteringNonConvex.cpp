#include "ScatteringNonConvex.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>
#include <iostream>

#include "Incidence.h"
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

void ScatteringNonConvex::ScatterLight(std::vector<Beam> &scaterredBeams)
{
	SplitLightToBeams();
	ScatterBeams(scaterredBeams);
}

void ScatteringNonConvex::PushBeamsToTree(Facet *facet, const PolygonArray &polygons)
{
	auto id = RecomputeTrackId(0, facet->index);
	m_splitting.inBeam.SetTracingParams(facet, 0, Location::In);
	m_splitting.outBeam.SetTracingParams(facet, 0, Location::Out);
	m_splitting.inBeam.id = id;
	m_splitting.outBeam.id = id;

	for (int j = 0; j < polygons.size; ++j)
	{
		Beam in = m_splitting.inBeam;
		Beam out = m_splitting.outBeam;
		const Polygon &pol = polygons.arr[j];
		in = pol;
		out = pol;
#ifdef _DEBUG // DEB
		in.pols.push_back(pol);
		out.pols.push_back(pol);
#endif
		Point3f p = pol.Center();
		in.opticalPath = 0;
		out.opticalPath = 0;
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
#ifdef _DEBUG // DEB
		in.ops.push_back(path);
		out.ops.push_back(path);
#endif
		in.AddOpticalPath(path);
		out.AddOpticalPath(path);
#ifdef _DEBUG // DEB
		in.dirs.push_back(in.direction);
		out.dirs.push_back(out.direction);
#endif
		m_propagatingBeams[m_treeSize++] = in;
		m_propagatingBeams[m_treeSize++] = out;
#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facet->in_normal, out);
#endif
	}
}

void ScatteringNonConvex::SplitByFacet(const Array<Facet*> &facets, int nCheckedFacets)
{

	PolygonArray resPolygons;
	IntersectWithFacet(facets, nCheckedFacets, resPolygons);

	if (resPolygons.size != 0)
	{
		Facet *facet = facets.elems[nCheckedFacets];

		m_splitting.SetBeams(Polygon());
		SetIncidentBeamOpticalParams(facet);
		PushBeamsToTree(facet, resPolygons);
	}
}

void ScatteringNonConvex::SplitLightToBeams()
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
#endif
	m_treeSize = 0;

	Array<Facet*> facets;
	SelectVisibleFacetsForLight(facets);

	for (int i = 0; i < facets.nElems; ++i)
	{
		SplitByFacet(facets, i);
	}
}

void ScatteringNonConvex::SelectVisibleFacetsForLight(Array<Facet *> &facetIDs)
{
	FindVisibleFacetsForLight(facetIDs);
	SortFacets_faster(m_incidentDir, facetIDs);
}

void ScatteringNonConvex::IntersectWithFacet(const Array<Facet*> &facets,
											 int nCheckedFacets,
											 PolygonArray &resFacets)
{
	Facet *facet = facets.elems[nCheckedFacets];

	if (facet->isOverlayedOut)
	{
		resFacets.Push(*facet);
	}
	else // facet is probably shadowed by others
	{
		const Point3f &normal = facet->ex_normal;
		CutPolygonByFacets(*facet, facets, nCheckedFacets, normal, normal,
						   m_incidentDir, resFacets);
	}
}

void ScatteringNonConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	FindVisibleFacets(beam, facets);

	Point3f dir = beam.direction;
	dir.d_param = beam.facet->in_normal.d_param;
	SortFacets_faster(dir, facets);
}

void ScatteringNonConvex::CutPolygonByFacets(const Polygon &pol,
											 const Array<Facet*> &facets, int size,
											 const Vector3f &polNormal,
											 const Vector3f &clipNormal,
											 const Vector3f &dir,
											 PolygonArray &pols)
{
	pols.Push(pol);

	// cut facet projections out of polygon one by one
	for (int i = 0; i < size; ++i)
	{
		m_polygonBuffer.Clear();

		while (pols.size != 0)
		{
			const Polygon &subj = pols.Pop();
			const Polygon &clip = *facets.elems[i];

			/// REF: объединить 2 первых аргумента и 2 вторых
			Difference(subj, polNormal, clip, clipNormal, dir, m_polygonBuffer);
		}

		if (m_polygonBuffer.size != 0)
		{
			for (int i = 0; i < m_polygonBuffer.size; ++i)
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
	const Point3f &n1 = beam.facet->ex_normal;
	const Point3f &n2 = beam.facet->in_normal;

	Array<Facet*> facets;
	SelectVisibleFacets(beam, facets);

	PolygonArray resultBeams;
	CutPolygonByFacets(beam, facets, facets.nElems, n1, n2,
					   -beam.direction, resultBeams);

	Beam tmp = beam;
	double path = m_splitting.ComputeOutgoingOpticalPath(tmp); // добираем оптический путь
	tmp.opticalPath += path;
#ifdef _DEBUG // DEB
	tmp.ops.push_back(path);
#endif
	for (int i = 0; i < resultBeams.size; ++i)
	{
		tmp.SetPolygon(resultBeams.arr[i]);
		scaterredBeams.push_back(tmp);
	}
}

void ScatteringNonConvex::SortFacets_faster(const Point3f &beamDir,
											Array<Facet*> &facets)
{
	if (facets.nElems == 0)
	{
		return;
	}

	int vertices[MAX_VERTEX_NUM];

	for (int i = 0; i < facets.nElems; ++i)
	{
		vertices[i] = FindClosestVertex(*facets.elems[i], beamDir);
	}

	int left = 0;
	int rigth = facets.nElems - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

	while (true)
	{
		int id = (left + rigth)/2;
		int i = left;
		int j = rigth;

		Point3f base = facets.elems[id]->arr[vertices[id]];

		while (i <= j)
		{
			Point3f vecB;
			double cosVN;

			do
			{
				vecB = base - facets.elems[i]->arr[vertices[i]];
				cosVN = DotProduct(vecB, beamDir);
				++i;
			}
			while (cosVN > FLT_EPSILON);
			--i;

			do
			{
				vecB = base - facets.elems[j]->arr[vertices[j]];
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

				Facet *temp_v = facets.elems[i];
				facets.elems[i] = facets.elems[j];
				facets.elems[j] = temp_v;

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

	for (int i = 1; i < facet.nVertices; ++i)
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

void ScatteringNonConvex::CutBeamByFacet(Facet *facet, Beam &beam,
										 PolygonArray &result)
{
	const Location &loc = beam.location;
	Facet *beamFacet = beam.facet;

	if (loc == Location::In && beamFacet->isOverlayedIn)
	{
		return;
	}

	const Point3f &facetNormal = (loc == Location::Out) ? -beamFacet->normal[loc]
														:  beamFacet->normal[loc];
	Difference(beam, beamFacet->normal[loc],
			   *facet, facetNormal, -beam.direction, result);

	if (result.size == 0) // beam is totaly swallowed by facet
	{
		beam.nVertices = 0;
	}
}

bool ScatteringNonConvex::IsOutgoingBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Out
			&& incidentBeam.nVertices != 0); // OPT: replace each other
}

int ScatteringNonConvex::FindFacetId(int facetId, const Array<Facet*> &arr)
{
	int i = 0;

	while ((facetId == arr.elems[i]->index) && (i < arr.nElems))
	{
		++i;
	}

	if (i == arr.nElems)
	{
		i = -1;
	}

	return i;
}

void ScatteringNonConvex::ScatterBeams(std::vector<Beam> &scaterredBeams)
{
#ifdef _DEBUG // DEB
	int count = 0;
#endif
	while (m_treeSize != 0)
	{
#ifdef _DEBUG // DEB
	++count;
#endif
		Beam beam = m_propagatingBeams[--m_treeSize];

		if (IsTerminalAct(beam) && beam.location == Location::Out) // REF, OPT: перенести проверку во все места, где пучок закидывается в дерево, чтобы пучки заранее не закидывались в него
		{
			CutExternalBeam(beam, scaterredBeams);
		}
		else
		{
			Array<Facet*> facets;
			SelectVisibleFacets(beam, facets);

			SplitBeamByVisibleFacets(facets, beam);

			if (IsOutgoingBeam(beam))
			{	// посылаем обрезанный всеми гранями внешний пучок на сферу
				double path = m_splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
				beam.opticalPath += path;
#ifdef _DEBUG // DEB
				beam.ops.push_back(path);
#endif
				scaterredBeams.push_back(beam);
			}
		}
	}
}

void ScatteringNonConvex::FindVisibleFacetsForLight(Array<Facet*> &facets)
{
	for (int i = 0; i < m_particle->nElems; ++i)
	{
		Facet *facet = m_particle->GetActualFacet(i);
		double cosA = DotProduct(m_incidentDir, facet->in_normal);

		if (cosA >= FLT_EPSILON) // beam incidents to this facet
		{
			facets.Add(facet);
		}
	}
}

bool ScatteringNonConvex::IsVisibleFacet(Facet *facet, const Beam &beam)
{
//	int loc = !beam.location;
	const Point3f &beamNormal = -beam.facet->normal[!beam.location];

	const Point3f &facetCenter = facet->center;
	const Point3f &beamCenter = beam.facet->center;
	Point3f vectorFromBeamToFacet = facetCenter - beamCenter;

	double cosBF = DotProduct(beamNormal, vectorFromBeamToFacet);
	return (cosBF >= EPS_ORTO_FACET);
}

void ScatteringNonConvex::FindVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	int begin = 0;
	int end = m_particle->nElems;

	if (m_particle->isAggregated && beam.location == Location::In)
	{
		m_particle->GetParticalFacetIdRange(beam.facet, begin, end);
	}

	for (int i = begin; i < end; ++i)
	{
		Facet *facet = m_particle->GetActualFacet(i);
		const Point3f &facetNormal = facet->normal[!beam.location];
		double cosFB = DotProduct(beam.direction, facetNormal);

		if (cosFB >= FLT_EPSILON) // beam incidents to this facet
		{
			if (IsVisibleFacet(facet, beam))
			{	// facet is in front of begin of beam
				facets.Add(facet);
			}
		}
	}
}

/// OPT: поменять все int и пр. параметры функций на ссылочные

/** TODO: придумать более надёжную сортировку по близости
 * (как вариант определять, что одна грань затеняют другую по мин. и макс.
 * удалённым вершинам, типа: "//" )
*/
void ScatteringNonConvex::SortFacets(const Point3f &beamDir, Array<Facet*> &facets)
{
	float distances[MAX_VERTEX_NUM];

	for (int i = 0; i < facets.nElems; ++i)
	{
		Facet *facet = facets.elems[i];
		distances[i] = CalcMinDistanceToFacet(facet, beamDir);
	}

	// sorting
	int left = 0;
	int rigth = facets.nElems - 1;

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

				Facet *temp_v = facets.elems[i];
				facets.elems[i] = facets.elems[j];
				facets.elems[j] = temp_v;

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

double ScatteringNonConvex::CalcMinDistanceToFacet(Polygon *facet,
												   const Point3f &beamDir)
{
	double dist = FLT_MAX;
	const Point3f *pol = facet->arr;
	Point3f point;
	Point3f dir = -beamDir;
	double dp = DotProduct(dir, beamDir);

	for (int i = 0; i < facet->nVertices; ++i)
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

	for (int i = 0; i < parts.size; ++i)
	{
		tmp = parts.arr[i];
		assert(m_treeSize < MAX_BEAM_REFL_NUM);
		m_propagatingBeams[m_treeSize++] = tmp;
	}
}

bool ScatteringNonConvex::FindRestOfBeamShape(Facet *facet, Beam &beam)
{
	m_polygonBuffer.Clear();
	CutBeamByFacet(facet, beam, m_polygonBuffer);

	bool isDivided = m_polygonBuffer.size > CLIP_RESULT_SINGLE;

	if (isDivided)
	{	// beam had divided by facet
		PushBeamPartsToTree(beam, m_polygonBuffer);
		beam.nVertices = 0;
	}
	else if (m_polygonBuffer.size == CLIP_RESULT_SINGLE)
	{
		beam = m_polygonBuffer.arr[0];
	}

	return isDivided;
}

bool ScatteringNonConvex::SplitBeamByVisibleFacets(const Array<Facet*> &facets, Beam &beam)
{
	bool isDivided = false;

	for (int i = 0; (i < facets.nElems) && !isDivided; ++i)// OPT: move this loop to SplitBeamByFacet
	{
		Facet *facet = facets.elems[i];

		Polygon beamShape;
		bool isIntersected = IncidentBeamToFacet(facet, beam, beamShape);

		if (isIntersected)
		{
			m_splitting.ComputeCosA(facet->ex_normal, beam.direction);
			m_splitting.SetBeams(beamShape);
#ifdef _DEBUG // DEB
			inBeam.pols = beam.pols;
			inBeam.pols.push_back(beamShape);
			outBeam.pols = beam.pols;
			outBeam.pols.push_back(beamShape);
#endif
			bool hasOutBeam = SetOpticalBeamParams(facet, beam);
			auto newId = RecomputeTrackId(beam.id, facet->index);

			if (hasOutBeam)
			{
				PushBeamToTree(m_splitting.outBeam, beam, newId, facet, Location::Out);
			}

			PushBeamToTree(m_splitting.inBeam, beam, newId, facet, Location::In);
			isDivided = FindRestOfBeamShape(facet, beam);
		}
	}
}

void ScatteringNonConvex::PushBeamsToBuffer(Facet *facet, const Beam &beam, bool hasOutBeam,
									   Beam &inBeam, Beam &outBeam,
									   std::vector<Beam> &passed)
{
	inBeam.id = beam.id;

	if (hasOutBeam)
	{
		outBeam.id = beam.id;
		outBeam.SetTracingParams(facet, beam.nActs+1, Location::Out);
		outBeam.id = RecomputeTrackId(outBeam.id, outBeam.facet->index);
		passed.push_back(outBeam);
	}

	inBeam.SetTracingParams(facet, beam.nActs+1, Location::In);
	inBeam.id = RecomputeTrackId(outBeam.id, outBeam.facet->index);
	passed.push_back(inBeam);
}

void ScatteringNonConvex::ScatterLight(const std::vector<std::vector<int>> &tracks,
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

//		for (int i = 1; i < track.size(); ++i)
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

//						Facet &facet = m_particle->GetActualFacet(facetID];
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
//}

//void ScatteringNonConvex::TraceFirstBeamFixedFacet(int facetID, bool &isIncident)
//{
//	isIncident = false;

//	IntArray facetIDs;
//	SelectVisibleFacetsForLight(facetIDs);

//	int index = FindFacetId(facetID, facetIDs);

//	for (int i = 0; (facetIDs.arr[i] < index) || (i < facetIDs.size); ++i)
//	{
//		int id = facetIDs.arr[i];
//		CutBeamByFacet();
//	}
//	if (index != -1)
//	{
//		SplitByFacet(facetIDs, index);
//		isIncident = true;
//	}
}
