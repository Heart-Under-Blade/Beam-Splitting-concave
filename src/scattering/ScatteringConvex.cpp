#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, const Light &incidentLight,
								   int maxActNo)
	: Scattering(particle, incidentLight, maxActNo)
{
}

void ScatteringConvex::PushBeamsToBuffer(Facet *facet, Splitting &splitting,
										 std::vector<Beam> &scatteredBeams)
{
	Track tr = m_originBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(0, facet->index);

	auto &beams = splitting.beams;
	beams.external.CopyTrack(tr);
	beams.external.SetLocation(false);

#ifdef MODE_FIXED_OR
	beams.external.dirs.push_back(beams.external.direction);
	beams.external.pols.push_back(beams.external);
#endif
	scatteredBeams.push_back(beams.external);

	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);

	PushBeamToTree(beams.internal);

#ifdef _CHECK_ENERGY_BALANCE
	ComputeFacetEnergy(facet->in_normal, m_splitting.beams.external);
#endif
}

void ScatteringConvex::SplitOriginBeam(std::vector<Beam> &scatteredBeams)
{
	m_visibleFacets.nElems = 0;
	SelectOriginVisibleFacets(m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		m_splitting.beams.SetBeams(*facet);
		ComputeOpticalBeamParams(facet, m_originBeam);
		PushBeamsToBuffer(facet, m_splitting, scatteredBeams);
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	FindVisibleFacets(beam, m_lightChecker, 0, m_particle->nElems, facets);
}
