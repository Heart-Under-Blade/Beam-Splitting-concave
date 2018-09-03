#include "TracingConcave.h"

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

TracingConcave::TracingConcave(Particle *particle, Light *incidentLight,
							   bool isOpticalPath, int interReflectionNumber)
	: Tracing(particle, incidentLight, isOpticalPath, interReflectionNumber)
{
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
		auto inID = inBeam.id;
		auto outID = outBeam.id;

		// set geometry of beam
		 inBeam.SetPolygon(polygons.arr[j]);
		outBeam.SetPolygon(polygons.arr[j]);

		// OPT: дублирует одиночный вызов в пред. ф-ции
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
	PolygonArray resPolygons;
	IntersectWithFacet(facetIDs, facetIndex, resPolygons);

#ifdef _DEBUG // DEB
	if (facetIDs.arr[facetIndex] == 74)
		int ff = 0;
#endif
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

	Point3f dir = beam.light.direction;
	dir.d_param = m_facets[beam.lastFacetID].in_normal.d_param;
	SortFacets_faster(dir, facetIDs);
}

void TracingConcave::CatchExternalBeam(const Beam &beam, std::vector<Beam> &scatteredBeams)
{
	const Point3f &normal = m_facets[beam.lastFacetID].ex_normal;
	const Point3f &normal1 = m_facets[beam.lastFacetID].in_normal;

	IntArray facetIds;
	SelectVisibleFacets(beam, facetIds);

	Polygon resultBeams[MAX_POLYGON_NUM];
	int resSize = 0;
	resultBeams[resSize++] = beam;

	Polygon diffFacets[MAX_POLYGON_NUM];

	// cut facet projections out of beam one by one
	for (int i = 0; i < facetIds.size; ++i)
	{
		int id = facetIds.arr[i];

		int diffSize = 0;

		while (resSize != 0)
		{
			Difference(resultBeams[--resSize], normal, m_facets[id],
					normal1, -beam.light.direction, diffFacets, diffSize);
		}

		if (diffSize == 0) // beam is totaly swallowed by facet
		{
			break;
		}

		for (int i = 0; i < diffSize; ++i)
		{
			resultBeams[resSize++] = diffFacets[i];
		}
	}

	Beam tmp = beam;
	tmp.opticalPath += fabs(FAR_ZONE_DISTANCE + tmp.D); // добираем оптический путь

	for (int i = 0; i < resSize; ++i)
	{
		tmp.SetPolygon(resultBeams[i]);
		scatteredBeams.push_back(tmp);
	}
}

void TracingConcave::SortFacets_faster(const Point3f &beamDir, IntArray &facetIDs)
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
			while (cosVN > EPS_COS_90);
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

int TracingConcave::FindClosestVertex(const Polygon &facet, const Point3f &beamDir)
{
	int closest = 0;

	for (int i = 1; i < facet.size; ++i)
	{
		Point3f v = facet.arr[closest] - facet.arr[i];
		double cosVD = DotProduct(v, beamDir);

		if (cosVD > EPS_COS_90)
		{
			closest = i;
		}
	}

	return closest;
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

void TracingConcave::CutBeamByFacet(int facetID, Beam &beam, bool &isDivided,
									Polygon *resultBeams, int &resultSize)
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
	Difference(beam, beamFacet.normal[loc],
			   m_facets[facetID], facetNormal, -beam.light.direction,
			   resultBeams, resultSize);

	if (resultSize == 0) // beam is totaly swallowed by facet
	{
		beam.size = 0;
	}
	else if (resultSize > CLIP_RESULT_SINGLE)
	{	// beam had divided by facet
		isDivided = true;
	}
}

bool TracingConcave::isExternalNonEmptyBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Out
			&& incidentBeam.size != 0); // OPT: replace each other
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
#ifdef _TRACK_ALLOW
	inBeam.id = beam.id;
	inBeam.locations = beam.locations;
#endif

	if (hasOutBeam)
	{
#ifdef _TRACK_ALLOW
		outBeam.locations = beam.locations;
		outBeam.id = beam.id;
#endif
		PushBeamToTree(outBeam, facetID, beam.level+1, Location::Out);
	}

	PushBeamToTree(inBeam, facetID, beam.level+1, Location::In);
}

void TracingConcave::TraceSecondaryBeams(std::vector<Beam> &scaterredBeams)
{
#ifdef _DEBUG // DEB
	int count = 0;
#endif
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];
#ifdef _DEBUG // DEB
		++count;

		vector<int> tr;
		Tracks::RecoverTrack(beam, m_particle->facetNum, tr);
		for (int gd : tr)
		{
			trackMapFile << gd << ' ';
		}
//		trackMapFile << bigIntegerToString(beam.id) << endl;
		trackMapFile << beam.id << endl;

		if (count == 52741)
			int fg = 0;
#endif
#ifdef _DEBUG // DEB
	if (beam.lastFacetID == 74 && beam.level == 0)
		int ff = 0;
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

		for (int i = 0; i < facetIds.size; ++i)// OPT: move this loop to TraceSecondaryBeamByFacet
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
			beam.opticalPath += fabs(FAR_ZONE_DISTANCE + beam.D); // добираем оптический путь
			scaterredBeams.push_back(beam);
		}
	}
}

void TracingConcave::SetOpticalBeamParams(int facetID, const Beam &incidentBeam,
										  Beam &inBeam, Beam &outBeam,
										  bool &hasOutBeam)
{
	const Point3f &dir = incidentBeam.light.direction;
	const Point3f &normal = m_facets[facetID].ex_normal;

	hasOutBeam = true;
	double cosIN = DotProduct(dir, normal);

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
#ifdef _DEBUG // DEB
			inBeam.dirs = incidentBeam.dirs;
			outBeam.dirs = incidentBeam.dirs;
#endif
			double cosI = DotProduct(-normal, dir);

			SetSloppingBeamParams_initial(dir, cosI, facetID, inBeam, outBeam);

			if (m_isOpticalPath)
			{
				CalcOpticalPath(cosIN, incidentBeam, inBeam, outBeam);
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
	int begin = 0;
	int end = m_particle->facetNum;

	if (m_particle->isAggregated && beam.location == Location::In)
	{
		m_particle->GetAggPartFacetIDRange(beam.lastFacetID, begin, end);
	}

	for (int i = begin; i < end; ++i)
	{
		const Point3f &facetNormal = m_facets[i].normal[!beam.location];
		double cosFB = DotProduct(beam.light.direction, facetNormal);

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

		while (resFacets.size != 0)
		{
			const Polygon &clip = m_facets[id];
			const Polygon &subj = resFacets.arr[--resFacets.size];
			Difference(subj, normal, clip, normal, m_incidentDir,
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
	Point3f point;
	Point3f dir = -beamDir;
	double dp = DotProduct(dir, beamDir);

	for (int i = 0; i < facet.size; ++i)
	{
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

void TracingConcave::TraceSecondaryBeamByFacet(Beam &beam, int facetID,
											   bool &isDivided)
{
	isDivided = false;
	Polygon intersected;
	bool hasIntersection = Intersect(facetID, beam, intersected);
#ifdef _DEBUG // DEB
	if (beam.lastFacetID == 8 && facetID == 65)
		int ff =0 ;
#endif
	if (hasIntersection)
	{
		Beam inBeam, outBeam;
		inBeam.SetPolygon(intersected);
		outBeam.SetPolygon(intersected);

		bool hasOutBeam;
		SetOpticalBeamParams(facetID, beam, inBeam, outBeam, hasOutBeam);
		inBeam.locations = beam.locations;
		outBeam.locations = beam.locations;
		PushBeamsToTree(beam, facetID, hasOutBeam, inBeam, outBeam);

		Polygon resultBeams[MAX_VERTEX_NUM];
		int resultSize = 0;
		CutBeamByFacet(facetID, beam, isDivided, resultBeams, resultSize);

		if (isDivided)
		{
			Beam tmp = beam;// OPT: try to replace 'tmp' to 'beam'

			for (int i = 0; i < resultSize; ++i)
			{
				tmp = resultBeams[i];
#ifdef _DEBUG // DEB
				if (m_treeSize >= MAX_BEAM_REFL_NUM)
					int f =0;
#endif
				assert(m_treeSize < MAX_BEAM_REFL_NUM);
				m_beamTree[m_treeSize++] = tmp;
			}
			beam.size = 0;
		}
		else if (resultSize == 1)
		{
			beam = resultBeams[0];
		}
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

void TracingConcave::SplitBeamByParticle(double beta, double gamma,
										 const std::vector<std::vector<int>> &tracks,
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

		for (size_t i = 1; i < track.size(); ++i)
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
//						CutBeamByFacet(facetID, beam, isDivided);
					}
				}
			}

			if (buffer.empty())
			{
				break;
			}

			for (const Beam &b : buffer)
			{	// добавляем прошедшие пучки
#ifdef _DEBUG // DEB
				if (m_treeSize >= MAX_BEAM_REFL_NUM)
					int f =0;
#endif
				assert(m_treeSize < MAX_BEAM_REFL_NUM);
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
// REF: replace 'wavefront' to something certain in your paper
	SelectVisibleFacetsForWavefront(facetIDs);
	int index = FindFacetID(facetID, facetIDs);

	if (index != -1)
	{
		TraceByFacet(facetIDs, index);
		isIncident = true;
	}
}
