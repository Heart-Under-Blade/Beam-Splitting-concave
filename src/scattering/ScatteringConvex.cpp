#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(const complex &refractiveIndex,
								   int maxActNo, double minBeamEnergy,
								   double farFresnelZone)
	: Scattering(refractiveIndex, maxActNo, minBeamEnergy, farFresnelZone)
{
}

void ScatteringConvex::SetWorkFacetsFromTracks()
{
	m_workFacets.nElems = 0;

	for (auto n : m_trackTreeNode->children)
	{
		auto facet = m_particle->GetActualFacet(n->m_facet->index);
		m_workFacets.Add(facet);
	}
}

void ScatteringConvex::Scatter(TrackNode *trackTree,
									std::vector<Beam> &scatteredBeams)
{
	m_scatteredBeams = &scatteredBeams;
	m_hasTracks = true;
	m_trackTreeNode = trackTree;
	SetWorkFacetsFromTracks();

	SplitStartBeam();
	// TODO: to complete
}

double ScatteringConvex::MeasureOpticalPath(const BeamInfo &info,
											const Point3f sourcePoint)
{
	double path = 0;
	Point3f dir = -info.beam->direction; // back direction
	bool loc = false;

	Point3f p1 = sourcePoint;
	Point3f p2;

	auto &tr = info.track;

	// back tracing
	for (int i = tr.size()-1; i > 0; --i)
	{
		Point3f &exNormal = m_particle->elems[tr[i]].actual.ex_normal;
		dir = m_splitting->ChangeBeamDirectionConvex(dir, exNormal, loc);

		Point3f &inNormal = m_particle->elems[tr[i-1]].actual.in_normal;
		p2 = Geometry::ProjectPointToPlane(p1, dir, inNormal);
		double len = Point3f::Length(p2 - p1);

		path += len;
		p1 = p2;
		loc = true;
	}

	lastPoint = p2;

	path *= real(m_splitting->GetRi());
	return path;
}

double ScatteringConvex::MeasureFullOpticalPath(const BeamInfo &info,
												const Point3f sourcePoint)
{
	double path = MeasureOpticalPath(info, sourcePoint);

	Point3f nFar1 = m_startBeam.direction;
	Point3f nFar2 = -info.beam->direction;
	path += m_farFresnelZone + Point3d::DotProduct(lastPoint, nFar1);
	path += fabs(Point3d::DotProduct(sourcePoint, nFar2) + m_farFresnelZone);

	return path;
}

void ScatteringConvex::SplitStartBeam()
{
	m_incidentEnergy = 0;

	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_startBeam, m_lightChecker, m_workFacets,
					  &m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];
		ComputeOpticalBeamParams(facet, m_startBeam, facet->ReversePolygon());
		ResolveBeams(m_startBeam, facet);
#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facet->in_normal, m_splitting->beams.external);
#endif
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam,
										   Array<Facet*> *facets)
{
	if (m_hasTracks)
	{
		SetWorkFacetsFromTracks();
	}

	FindVisibleFacets(beam, m_lightChecker, m_workFacets, facets);
}

void ScatteringConvex::SplitSecondaryBeams()
{
	m_isBeamInside = true;
#ifdef _DEBUG // DEB
	int cc = 0;
#endif
//	bool isEmpty = false;

	while (m_secondaryBeams->size != 0)
	{
#ifdef _DEBUG // DEB
		/*cout <<*/ ++cc /*<< endl*/;
		if (m_scatteredBeams->size() == 15)
			int ddd = 0;
#endif
		Beam beam = m_secondaryBeams->Pop();
		SplitBeamByVisibleFacets(beam);
	}
}

void ScatteringConvex::ResolveExternalBeam(const Beam &beam)
{
	m_scatteredBeams->push_back(beam);
}
