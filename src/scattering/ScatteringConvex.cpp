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

	splitting.external.CopyTrack(tr);
	splitting.external.SetLocation(false);
	scatteredBeams.push_back(splitting.external);

	splitting.internal.CopyTrack(tr);
	splitting.internal.SetLocation(true);

	if (IsTerminalAct(splitting.internal))
	{
		if (!splitting.internal.isInside)
		{
			ReleaseBeam(splitting.internal);
		}
	}
	else
	{
		PushBeamToTree(splitting.internal);
	}

#ifdef _CHECK_ENERGY_BALANCE
	ComputeFacetEnergy(facet->in_normal, m_splitting.external);
#endif
}

void ScatteringConvex::SplitOriginBeam(std::vector<Beam> &scatteredBeams)
{
	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originBeam, m_lightChecker, 0, m_particle->nElems, m_visibleFacets);

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
