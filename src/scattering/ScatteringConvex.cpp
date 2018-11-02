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

	splitting.outBeam.CopyTrack(tr);
	splitting.outBeam.SetLocation(false);
	scatteredBeams.push_back(splitting.outBeam);

	splitting.inBeam.CopyTrack(tr);
	splitting.inBeam.SetLocation(true);

	if (IsTerminalAct(splitting.inBeam))
	{
		if (!splitting.inBeam.isInside)
		{
			ReleaseBeam(splitting.inBeam);
		}
	}
	else
	{
		PushBeamToTree(splitting.inBeam);
	}

#ifdef _CHECK_ENERGY_BALANCE
	ComputeFacetEnergy(facet->in_normal, m_splitting.outBeam);
#endif
}

void ScatteringConvex::SplitOriginBeam(std::vector<Beam> &scatteredBeams)
{
	m_visibleFacets.nElems = 0;
	SelectOriginVisibleFacets(m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		m_splitting.SetBeams(*facet);
		ComputeOpticalBeamParams(facet, m_originBeam);
		PushBeamsToBuffer(facet, m_splitting, scatteredBeams);
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	FindVisibleFacets(beam, m_lightChecker, 0, m_particle->nElems, facets);
}
