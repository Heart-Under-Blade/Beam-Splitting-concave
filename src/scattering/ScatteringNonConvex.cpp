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
										 const Light &incidentLight,
										 int maxActNo,
										 const complex &refractiveIndex)
	: Scattering(particle, incidentLight, maxActNo, refractiveIndex)
{
}

ScatteringNonConvex::~ScatteringNonConvex()
{
}

void ScatteringNonConvex::ScatterLight(TrackNode *trackTree,
									   std::vector<Beam> &scatteredBeams)
{
	m_hasTracks = true;
	m_trackTreeNode = trackTree;
	SplitOriginalBeam(scatteredBeams);
	// TODO: to complete
}

void ScatteringNonConvex::SplitOriginalBeam(std::vector<Beam> &scatteredBeams)
{
	m_incidentEnergy = 0;

	auto &allFacets = m_workFacets;

	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, allFacets,
					  m_visibleFacets);
	SortFacetsByDistance(m_originalBeam.direction, m_visibleFacets);

	auto &visiblePart = m_intersectionBuffer;

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		if (!m_hasTracks || m_trackTreeNode->FindNode(facet->index) != nullptr)
		{
			visiblePart.Clear();
			bool isFound = FindVisiblePartOfFacet(m_visibleFacets, i, visiblePart);

			if (isFound)
			{
				ComputeOpticalBeamParams(facet, m_originalBeam, *facet);
				PushBeamsToTree(facet, m_splitting.beams, visiblePart,
								scatteredBeams);
			}
		}
	}
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
	if (Scattering::IsFinalAct(beam))
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
											  Array<Facet*> &visibleFacets)
{
	if (m_particle->isAggregated && beam.isInside)
	{
		m_workFacets.nElems = 0;
		m_particle->GetPartByFacet(beam.facet, m_workFacets);
	}

	FindVisibleFacets(beam, m_beamChecker, m_workFacets, visibleFacets);

	Point3f dir = beam.direction;
	dir.d_param = beam.facet->in_normal.d_param;
	SortFacetsByDistance(dir, visibleFacets);
}

bool ScatteringNonConvex::FindVisiblePartOfFacet(const Array<Facet*> &facets,
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
	if (Scattering::IsFinalAct(beam))
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
//		tmp.ops.push_back(path);
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

bool ScatteringNonConvex::FindRestOfBeamShape(Facet *facet, const Beam &beam,
											  PolygonStack &rest)
{
	bool isSwallowed = false;

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
void ScatteringNonConvex::PushBeamPartsToBuffer(const Beam &beam,
												const PolygonStack &parts)
{
	Beam tmp = beam; // OPT: try to replace 'tmp' to 'beam'

	for (int i = 0; i < parts.nPolygons; ++i)
	{
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

bool ScatteringNonConvex::IsFinalAct(const Beam &beam)
{
	return Scattering::IsFinalAct(beam) || (!beam.isInside && beam.nVertices != 0);
}
