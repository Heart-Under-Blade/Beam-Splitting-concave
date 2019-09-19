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

ScatteringNonConvex::ScatteringNonConvex(const complex &refractiveIndex,
										 int maxActNo, double minBeamEnergy,
										 double farFresnelZone)
	: Scattering(refractiveIndex, maxActNo, minBeamEnergy, farFresnelZone)
{
}

ScatteringNonConvex::~ScatteringNonConvex()
{
}

void ScatteringNonConvex::Scatter(TrackNode *trackTree,
								  std::vector<Beam> &scatteredBeams)
{
	m_scatteredBeams = &scatteredBeams;
	m_hasTracks = true;
	m_trackTreeNode = trackTree;
	SplitStartBeam();
	// TODO: to complete
}

double ScatteringNonConvex::MeasureOpticalPath(const BeamInfo &info,
											   const Point3f sourcePoint)
{
	Beam *beam = info.beam;
	double path = 0;
	Point3f dir = -beam->direction; // back direction
	bool loc = false;
	bool nextLoc;

	Point3f p1 = sourcePoint;
	Point3f p2;

	auto &tr = info.track;

	// back tracing
	for (int i = tr.size()-1; i > 0; --i)
	{
		nextLoc = beam->IsInsideAtAct(i-1);

		Point3f &exNormal = m_particle->elems[tr[i]].actual.ex_normal;
		dir = m_splitting->ChangeBeamDirection(dir, exNormal, loc, nextLoc);

		Point3f &inNormal = m_particle->elems[tr[i-1]].actual.in_normal;
		p2 = Geometry::ProjectPointToPlane(p1, dir, inNormal);
		double len = Point3f::Length(p2 - p1);

		if (nextLoc)
		{	// add internal path only
#ifdef _TEST
			len *= sqrt(real(m_splitting->GetRi()));
#else
			double cosA = Point3f::DotProduct(dir, exNormal);
			double reRi = m_splitting->ComputeEffectiveReRi(cosA);
			len *= sqrt(reRi);
#endif
		}
#ifdef _DEBUG // DEB
//		lens.push_back(len);
#endif
		path += len;
		p1 = p2;
		loc = nextLoc;
	}

#ifdef _DEBUG // DEB
//	path *= real(m_splitting->GetRi());
//	Point3f nFar1 = m_incidentDir;
//	Point3f nFar2 = -beam.direction;
//	double dd1 = m_splitting->FAR_ZONE_DISTANCE + DotProductD(p2, nFar1);
//	double dd2 = fabs(DotProductD(sourcePoint, nFar2) + m_splitting->FAR_ZONE_DISTANCE);
//	path += dd1;
//	path += dd2;
//	if (fabs(path - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}

void ScatteringNonConvex::SplitStartBeam()
{
	m_incidentEnergy = 0;

	auto &allFacets = m_workFacets;

	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_startBeam, m_lightChecker, allFacets,
					  &m_visibleFacets);
	SortFacetsByDistance(m_startBeam.direction, &m_visibleFacets);

	auto &visiblePart = m_splittedBeams;

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		if (!m_hasTracks || m_trackTreeNode->FindNode(facet->index) != nullptr)
		{
			visiblePart.Clear();
			bool isFound = FindVisibleFacetPart(m_visibleFacets, i, visiblePart);

			if (isFound)
			{
				ComputeOpticalBeamParams(facet, m_startBeam, *facet);
				ResolveBeamParts(facet, visiblePart);
			}
		}
	}
}

void ScatteringNonConvex::ResolveBeamParts(Facet *facet,
										   const PolygonStack &parts)
{
	auto &beams = m_splitting->beams;
	Track beamTrack;
	beamTrack.Update(facet);
	beamTrack.RefreshId(facet->index, m_particle->nElems);

	if (!IsFinalAct(beams.internal))
	{
		beams.internal.CopyTrack(beamTrack);
		m_secondaryBeams = &m_internalBeams;
		beams.internal.SetIsInside(true);
		StackBeamParts(beams.internal, parts);
	}

	beams.external.CopyTrack(beamTrack);

	if (!IsFinalAct(beams.external))
	{
		m_secondaryBeams = &m_externalBeams;
		beams.internal.SetIsInside(false);
		StackBeamParts(beams.external, parts);
	}
	else
	{
		Beam beam = beams.external;

		for (int i = 0; i < parts.nPolygons; ++i)
		{
			beam.SetPolygon(parts.polygons[i]);
			SplitExternalBeam(beam);
			ReleaseFinalBeams(beam);
		}
	}

#ifdef _CHECK_ENERGY_BALANCE
	for (int j = 0; j < parts.nPolygons; ++j)
	{
		ComputeFacetEnergy(facet->in_normal, parts.polygons[j]);
	}
#endif
}

void ScatteringNonConvex::ResolveBeams(Beam &parentBeam, Facet *facet)
{
	Scattering::ResolveBeams(parentBeam, facet);

	m_isDivided = false;

	m_beamParts.Clear();
	PolygonResult pr = FindRestOfBeamShape(facet, parentBeam, m_beamParts);

	if (pr == PolygonResult::Multiple)
	{
		m_isDivided = true;
		StackBeamParts(parentBeam, m_beamParts);
		parentBeam.Clear();
	}
	else if (pr == PolygonResult::Single)
	{
		parentBeam = m_beamParts.polygons[0];

		if (m_isBeamInside)
		{
			StackBeam(parentBeam);
		}
		else
		{
			ResolveExternalBeam(parentBeam);
		}
	}
	else if (pr == PolygonResult::None) // beam is totally overlayed by facet
	{
		parentBeam.Clear();
	}
}

void ScatteringNonConvex::SelectVisibleFacets(const Beam &beam,
											  Array<Facet *> *visibleFacets)
{
	if (m_particle->isAggregated && m_isBeamInside)
	{
		m_workFacets.nElems = 0;
		m_particle->GetPartByFacet(beam.facet, m_workFacets);
	}

	FindVisibleFacets(beam, m_beamChecker, m_workFacets, visibleFacets);

	Point3f dir = beam.direction;
	dir.d_param = beam.facet->in_normal.d_param;
	SortFacetsByDistance(dir, visibleFacets);
}

bool ScatteringNonConvex::FindVisibleFacetPart(const Array<Facet*> &facets,
											   int facetIndex,
											   PolygonStack &pols)
{
	Facet *facet = facets.elems[facetIndex];

	if (facetIndex == 0 || !facet->isOverlayedOut)
	{
		pols.Push(*facet);
	}
	else // cut facet by shadows of others
	{
		const Point3f &normal = facet->ex_normal;
		CutPolygonByFacets(*facet, facets, facetIndex, normal, normal,
						   m_startBeam.direction, pols);
	}

	return pols.nPolygons != 0;
}

void ScatteringNonConvex::CutPolygonByFacets(const Polygon &pol,
											 const Array<Facet*> &facets,
											 int lastFacetIndex,
											 const Vector3f &polNormal,
											 const Vector3f &clipNormal,
											 const Vector3f &dir,
											 PolygonStack &pols)
{
	pols.Push(pol);
	m_beamParts.nPolygons = 1;

	// cut facet projections out of polygon one by one
	for (int i = 0; (i < lastFacetIndex && m_beamParts.nPolygons != 0); ++i)
	{
		m_beamParts.Clear();

		while (pols.nPolygons != 0)
		{
			const Polygon &subj = pols.Pop();
			const Polygon &clip = *facets.elems[i];

			/// REF: объединить 2 первых аргумента и 2 вторых
			Geometry::DifferPolygons(subj, polNormal, clip, clipNormal,
									 dir, m_beamParts);
		}

		for (int i = 0; i < m_beamParts.nPolygons; ++i)
		{
			pols.Push(m_beamParts.polygons[i]);
		}
	}
}

void ScatteringNonConvex::SplitExternalBeam(const Beam &beam)
{
	const Point3f &n1 = beam.facet->ex_normal;
	const Point3f &n2 = beam.facet->in_normal;

	m_visibleFacets.nElems = 0;
	SelectVisibleFacets(beam, &m_visibleFacets);

	m_splittedBeams.Clear();
	CutPolygonByFacets(beam, m_visibleFacets, m_visibleFacets.nElems, n1, n2,
					   -beam.direction, m_splittedBeams);
//	double path = m_splitting->ComputeOutgoingOpticalPath(tmp); // добираем оптический путь
//	tmp.opticalPath += path;
#ifdef _DEBUG // DEB
//	tmp.ops.push_back(path);
#endif
}

void ScatteringNonConvex::ReleaseFinalBeams(Beam beam)
{
	for (int i = 0; i < m_splittedBeams.nPolygons; ++i)
	{
#ifdef MODE_FIXED_OR
		tmp.pols.push_back(tmp);
#endif
		beam.SetPolygon(m_splittedBeams.polygons[i]);
		m_scatteredBeams->push_back(beam);
	}
}

void ScatteringNonConvex::ResolveExternalBeam(const Beam &beam)
{
	if (IsFinalAct(beam))
	{
		if (beam.nVertices != 0)
		{
			SplitExternalBeam(beam);
			ReleaseFinalBeams(beam);
		}
	}
	else
	{
		m_externalBeams.Push(beam);
	}
}

void ScatteringNonConvex::SortFacetsByDistance(const Vector3f &beamDir,
											   Array<Facet*> *facets) const
{
	if (facets->nElems == 0)
	{
		return;
	}

	int vertices[MAX_VERTEX_NUM];

	for (int i = 0; i < facets->nElems; ++i)
	{
		vertices[i] = FindClosestVertex(*facets->elems[i], beamDir);
	}

	int left = 0;
	int rigth = facets->nElems - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

	while (true)
	{
		int id = (left + rigth)/2;
		int i = left;
		int j = rigth;

		Point3f base = facets->elems[id]->vertices[vertices[id]];

		while (i <= j)
		{
			Point3f vecB;
			double cosVN;

			do
			{
				vecB = base - facets->elems[i]->vertices[vertices[i]];
				cosVN = Point3f::DotProduct(vecB, beamDir);
				++i;
			}
			while (cosVN > FLT_EPSILON);
			--i;

			do
			{
				vecB = base - facets->elems[j]->vertices[vertices[j]];
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

				Facet *temp_v = facets->elems[i];
				facets->elems[i] = facets->elems[j];
				facets->elems[j] = temp_v;

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

	for (int i = 1; i < facet.nVertices; ++i)
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

PolygonResult ScatteringNonConvex::FindRestOfBeamShape(Facet *facet,
													   const Beam &beam,
													   PolygonStack &rest)
{
	PolygonResult pr = PolygonResult::Uncut;

	if (!m_isBeamInside || beam.facet->isOverlayedIn)
	{
		Geometry::DifferPolygons(beam, beam.facet->normal[!m_isBeamInside],
				*facet, beam.facet->in_normal, -beam.direction, rest);

		if (rest.nPolygons == PolygonResult::None)
		{
			pr = PolygonResult::None;
		}
		else if (rest.nPolygons == PolygonResult::Single)
		{
			pr = PolygonResult::Single;
		}
		else if (rest.nPolygons > PolygonResult::Single)
		{
			pr = PolygonResult::Multiple;
		}
	}

	return pr;
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

	for (int i = 0; i < facet->nVertices; ++i)
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
void ScatteringNonConvex::StackBeamParts(Beam beam, const PolygonStack &parts)
{
	for (int i = 0; i < parts.nPolygons; ++i)
	{
		beam.SetPolygon(parts.polygons[i]);
		assert(m_secondaryBeams->size < MAX_BEAM_NUM);
#ifdef MODE_FIXED_OR
		if (!m_isDivided)
		{
			tmp.dirs.push_back(tmp.direction);
			tmp.pols.push_back(tmp);
		}
#endif
		m_secondaryBeams->Push(beam);
	}
}

bool ScatteringNonConvex::isFinalFacet(int index, Array<Facet*> &facets)
{
	if (m_isDivided)
	{
		m_isDivided = false;
		return true;
	}
	else
	{
		return Scattering::isFinalFacet(index, facets);
	}
}


void ScatteringNonConvex::SplitSecondaryBeams()
{
	m_isBeamInside = true;
	bool isEmpty = false;

	while (!isEmpty)
	{
		if (m_secondaryBeams->size != 0)
		{
			Beam beam = m_secondaryBeams->Pop();
			SplitBeamByVisibleFacets(beam);
		}
		else
		{
			if (m_internalBeams.size != 0)
			{
				m_secondaryBeams = &m_internalBeams;
				m_isBeamInside = true;
			}
			else if (m_externalBeams.size != 0)
			{
				m_secondaryBeams = &m_externalBeams;
				m_isBeamInside = false;
			}
			else
			{
				isEmpty = true;
			}
		}
	}
}
