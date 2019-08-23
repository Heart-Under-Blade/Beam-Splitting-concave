#include "ScatteringNonConvex.h"

#include "macro.h"
#include <tgmath.h>
#include <assert.h>
#include <iostream>

#include "BigIntegerLibrary.hh"

//#ifdef _DEBUG // DEB
//#include "Tracer.h"
//#endif

#ifdef _DEBUG
using namespace std;
ofstream trackMapFile("tracks_deb.dat", ios::out);
#endif

using namespace std;

ScatteringNonConvex::ScatteringNonConvex(Particle *particle,
										 const Light &incidentLight, int maxActNo)
	: Scattering(particle, incidentLight, maxActNo)
{
}

ScatteringNonConvex::~ScatteringNonConvex()
{
<<<<<<< HEAD
//	m_particle->Rotate(beta, gamma, 0);
	SplitLightToBeams();
	SplitBeams(scaterredBeams);
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}

void ScatteringNonConvex::SplitOriginalBeam(std::vector<Beam> &externalBeams)
{
	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, 0, m_particle->nElems,
					  m_visibleFacets);
	SortFacetsByDistance(m_originalBeam.direction, m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
<<<<<<< HEAD
		Beam in = inBeam;
		Beam out = outBeam;
		const Polygon &pol = polygons.arr[j];
		in = pol;
		out = pol;
#ifdef _DEBUG // DEB
//		in.pols.push_back(pol);
//		out.pols.push_back(pol);
#endif
		Point3f p = pol.Center();
		in.opticalPath = 0;
		out.opticalPath = 0;
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
#ifdef _DEBUG // DEB
//		in.ops.push_back(path);
//		out.ops.push_back(path);
#endif
		in.AddOpticalPath(path);
		out.AddOpticalPath(path);
#ifdef _DEBUG // DEB
//		in.dirs.push_back(in.direction);
//		out.dirs.push_back(out.direction);
=======
		m_intersectionBuffer.Clear();
		bool isIntersected = FindLightedFacetPolygon(m_visibleFacets, i,
													 m_intersectionBuffer);
		if (isIntersected)
		{
			Facet *facet = m_visibleFacets.elems[i];

			m_splitting.SetBeams(*facet);
			ComputeOpticalBeamParams(facet, m_originalBeam);
			PushBeamsToTree(facet, m_splitting.beams, m_intersectionBuffer,
							externalBeams);
		}
	}

#ifdef _DEBUG // DEB
//	Beam b = m_propagatingBeams[18];
//	m_propagatingBeams[18] = m_propagatingBeams[17];
//	m_propagatingBeams[17] = b;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#endif
}

void ScatteringNonConvex::PushBeamsToTree(Facet *facet, BeamPair<Beam> &beams,
										  const PolygonStack &polygons,
										  std::vector<Beam> &scatteredBeams)
{
	Track beamTrack;
	beamTrack.Update(facet);
	beamTrack.RecomputeTrackId(facet->index, m_particle->nElems);

	beams.internal.CopyTrack(beamTrack);
	beams.internal.SetLocation(true);
	PushBeamToBuffer(beams.internal, polygons, scatteredBeams);

	beams.external.CopyTrack(beamTrack);
	beams.external.SetLocation(false);
	PushBeamToBuffer(beams.external, polygons, scatteredBeams);

#ifdef _CHECK_ENERGY_BALANCE
	for (int j = 0; j < polygons.nPolygons; ++j)
	{
		ComputeFacetEnergy(facet->in_normal, polygons.polygons[j]);
	}
#endif
}

void ScatteringNonConvex::PushBeamToBuffer(Beam &beam, const PolygonStack &beamParts,
										   std::vector<Beam> &scatteredBeams)
{
	if (Scattering::IsTerminalAct(beam))
	{
		if (!beam.isInside)
		{
			ReleaseBeam(beam);
		}
	}
	else
	{
		PushBeamPartsToBuffer(beam, beamParts);
	}
}

void ScatteringNonConvex::PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
											bool hasOutBeam)
{
	Scattering::PushBeamsToBuffer(parentBeam, facet, hasOutBeam);

	m_differenceBuffer.Clear();
	bool isOverlayed = FindRestOfBeamShape(facet, parentBeam, m_differenceBuffer);

	m_isDivided = m_differenceBuffer.nPolygons > CLIP_RESULT_SINGLE;

	if (m_isDivided)
	{
		PushBeamPartsToBuffer(parentBeam, m_differenceBuffer);
		parentBeam.nVertices = 0;
	}
	else if (m_differenceBuffer.nPolygons == CLIP_RESULT_SINGLE)
	{
		parentBeam = m_differenceBuffer.polygons[0];
	}
	else if (isOverlayed) // beam is totally overlayed by facet
	{
		parentBeam.nVertices = 0;
	}
}

void ScatteringNonConvex::SelectVisibleFacets(const Beam &beam,
											  Array<Facet*> &facets)
{
	int begin = 0;
	int end = m_particle->nElems;

	if (m_particle->isAggregated && beam.isInside)
	{
		m_particle->GetParticalFacetIdRange(beam.facet, begin, end);
	}

	FindVisibleFacets(beam, m_beamChecker, begin, end, facets);

	Point3f dir = beam.direction;
	dir.d_param = beam.facet->in_normal.d_param;
	SortFacetsByDistance(dir, facets);
}

bool ScatteringNonConvex::FindLightedFacetPolygon(const Array<Facet*> &facets,
												  int nCheckedFacets,
												  PolygonStack &pols)
{
	Facet *facet = facets.elems[nCheckedFacets];

	if (nCheckedFacets == 0 || !facet->isOverlayedOut)
	{
		pols.Push(*facet);
	}
	else // cut facet by shadows of others
	{
		const Point3f &normal = facet->ex_normal;
		CutPolygonByFacets(*facet, facets, nCheckedFacets, normal, normal,
						   m_originalBeam.direction, pols);
	}

	return pols.nPolygons != 0;
}

void ScatteringNonConvex::CutPolygonByFacets(const Polygon &pol,
											 const Array<Facet*> &facets, int size,
											 const Vector3f &polNormal,
											 const Vector3f &clipNormal,
											 const Vector3f &dir,
											 PolygonStack &pols)
{
	pols.Push(pol);
	m_differenceBuffer.nPolygons = 1;

	// cut facet projections out of polygon one by one
	for (int i = 0; (i < size && m_differenceBuffer.nPolygons != 0); ++i)
	{
		m_differenceBuffer.Clear();

		while (pols.nPolygons != 0)
		{
			const Polygon &subj = pols.Pop();
			const Polygon &clip = *facets.elems[i];

			/// REF: объединить 2 первых аргумента и 2 вторых
			Geometry::DifferPolygons(subj, polNormal, clip, clipNormal,
									 dir, m_differenceBuffer);
		}

		for (int i = 0; i < m_differenceBuffer.nPolygons; ++i)
		{
			pols.Push(m_differenceBuffer.polygons[i]);
		}
	}
}

void ScatteringNonConvex::ReleaseBeam(Beam &beam)
{
	if (Scattering::IsTerminalAct(beam))
	{
		const Point3f &n1 = beam.facet->ex_normal;
		const Point3f &n2 = beam.facet->in_normal;

		m_visibleFacets.nElems = 0;
		SelectVisibleFacets(beam, m_visibleFacets);

		m_intersectionBuffer.Clear();
		CutPolygonByFacets(beam, m_visibleFacets, m_visibleFacets.nElems, n1, n2,
						   -beam.direction, m_intersectionBuffer);
		Beam tmp = beam;
//		double path = m_splitting.ComputeOutgoingOpticalPath(tmp); // добираем оптический путь
//		tmp.opticalPath += path;
#ifdef _DEBUG // DEB
<<<<<<< HEAD
//	tmp.ops.push_back(path);
=======
//		tmp.ops.push_back(path);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#endif
		for (int i = 0; i < m_intersectionBuffer.nPolygons; ++i)
		{
#ifdef MODE_FIXED_OR
			tmp.pols.push_back(tmp);
#endif
			tmp.SetPolygon(m_intersectionBuffer.polygons[i]);
			m_scatteredBeams->push_back(tmp);
		}
	}
	else
	{
		Scattering::ReleaseBeam(beam);
	}
}

void ScatteringNonConvex::SortFacetsByDistance(const Vector3f &beamDir,
											   Array<Facet*> &facets) const
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

		Point3f base = facets.elems[id]->vertices[vertices[id]];

		while (i <= j)
		{
			Point3f vecB;
			double cosVN;

			do
			{
				vecB = base - facets.elems[i]->vertices[vertices[i]];
				cosVN = Point3f::DotProduct(vecB, beamDir);
				++i;
			}
			while (cosVN > FLT_EPSILON);
			--i;

			do
			{
				vecB = base - facets.elems[j]->vertices[vertices[j]];
				cosVN = Point3f::DotProduct(vecB, beamDir);
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

int ScatteringNonConvex::FindClosestVertex(const Polygon &facet,
										   const Point3f &beamDir) const
{
	int closest = 0;

<<<<<<< HEAD
	for (unsigned i = 1; i < facet.nVertices; ++i)
=======
	for (int i = 1; i < facet.nVertices; ++i)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	{
		Point3f v = facet.vertices[closest] - facet.vertices[i];
		double cosVD = Point3f::DotProduct(v, beamDir);

		if (cosVD > FLT_EPSILON)
		{
			closest = i;
		}
	}

	return closest;
}

<<<<<<< HEAD
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
		beam.nVertices = 0;
	}
}

bool ScatteringNonConvex::IsOutgoingBeam(Beam &incidentBeam)
{
	return (incidentBeam.location == Location::Out
			&& incidentBeam.nVertices != 0); // OPT: replace each other
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
//	int count = 0;
#endif
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];

#ifdef _DEBUG // DEB
//	++count;
#endif
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

				if (intersection.nVertices >= MIN_VERTEX_NUM)
				{
					isDivided = SplitBeamByFacet(intersection, facetId, beam);
				}
			}

			if (IsOutgoingBeam(beam))
			{	// посылаем обрезанный всеми гранями внешний пучок на сферу
				double path = m_splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
				beam.opticalPath += path;
#ifdef _DEBUG // DEB
//	beam.ops.push_back(path);
#endif
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
#ifdef _DEBUG // DEB
//	inBeam.ops = incidentBeam.ops;
//	outBeam.ops = incidentBeam.ops;
//	inBeam.ops.push_back(path);
//	outBeam.ops.push_back(path);
#endif
				path += incidentBeam.opticalPath;
				inBeam.AddOpticalPath(path);
				outBeam.AddOpticalPath(path);
			}
		}
	}
=======
bool ScatteringNonConvex::FindRestOfBeamShape(Facet *facet, const Beam &beam,
											  PolygonStack &rest)
{
	bool isSwallowed = false;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	if (!beam.isInside || beam.facet->isOverlayedIn)
	{
		Geometry::DifferPolygons(beam, beam.facet->normal[!beam.isInside],
				*facet, beam.facet->in_normal, -beam.direction, rest);

		isSwallowed = (rest.nPolygons == 0);
	}

	return isSwallowed;
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
	const Point3f *pol = facet->vertices;
	Point3f point;
	Point3f dir = -beamDir;
	double dp = Point3f::DotProduct(dir, beamDir);

<<<<<<< HEAD
	for (int i = 0; i < facet.nVertices; ++i)
=======
	for (int i = 0; i < facet->nVertices; ++i)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	{
		/// REF: заменить на сущ. фуyкцию ProjectPointToPlane
		// measure dist
		double t = Point3f::DotProduct(pol[i], beamDir);
		t = t + beamDir.d_param;
		t = t/dp;
		point = pol[i] - (dir * t);
		double newDist = sqrt(Point3f::Norm(point - pol[i]));

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
void ScatteringNonConvex::PushBeamPartsToBuffer(const Beam &beam,
												const PolygonStack &parts)
{
	Beam tmp = beam; // OPT: try to replace 'tmp' to 'beam'

	for (int i = 0; i < parts.nPolygons; ++i)
	{
<<<<<<< HEAD
		tmp = parts.arr[i];
		assert(m_treeSize < MAX_BEAM_REFL_NUM);
		m_beamTree[m_treeSize++] = tmp;
	}
}

template<class T>
void ScatteringNonConvex::PushBeamToTree(Beam &beam, const Beam &oldBeam,
										 const T &newId, int facetId,
										 Location loc)
{
	beam.id = newId;
	beam.locations = oldBeam.locations;
#ifdef _DEBUG // DEB
//	beam.dirs = oldBeam.dirs;
#endif
	Scattering::PushBeamToTree(beam, facetId, oldBeam.nActs+1, loc);
}

bool ScatteringNonConvex::SplitBeamByFacet(const Polygon &intersection,
										   int facetId, Beam &beam)
{
	Facet &facet = m_facets[facetId];

	Beam inBeam, outBeam;
	inBeam.SetPolygon(intersection);
	outBeam.SetPolygon(intersection);
#ifdef _DEBUG // DEB
//	inBeam.pols = beam.pols;
//	outBeam.pols = beam.pols;
//	inBeam.pols.push_back(intersection);
//	outBeam.pols.push_back(intersection);
#endif

	bool hasOutBeam = SetOpticalBeamParams(facet, beam, inBeam, outBeam);

	auto newId = RecomputeTrackId(beam.id, facetId);

	PushBeamToTree(inBeam, beam, newId, facetId, Location::In);

	if (hasOutBeam)
=======
		tmp = parts.polygons[i];
		assert(m_treeSize < MAX_BEAM_NUM);
#ifdef MODE_FIXED_OR
		if (!m_isDivided)
		{
			tmp.dirs.push_back(tmp.direction);
			tmp.pols.push_back(tmp);
		}
#endif
		m_propagatingBeams[m_treeSize++] = tmp;
	}
}

bool ScatteringNonConvex::isTerminalFacet(int index, Array<Facet*> &facets)
{
	if (m_isDivided)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	{
		m_isDivided = false;
		return true;
	}
<<<<<<< HEAD

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
=======
	else
	{
		return Scattering::isTerminalFacet(index, facets);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	}
}

bool ScatteringNonConvex::IsTerminalAct(const Beam &beam)
{
	return Scattering::IsTerminalAct(beam) || (!beam.isInside && beam.nVertices != 0);
}
