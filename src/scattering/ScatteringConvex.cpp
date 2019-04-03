#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, const Light &incidentLight,
								   int maxActNo)
	: Scattering(particle, incidentLight, maxActNo)
{
}

void ScatteringConvex::PushBeamsToBuffer(Facet *facet, BeamPair<Beam> &beams,
										 std::vector<Beam> &scatteredBeams)
{
	Track tr = m_originalBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(0, facet->index);

	beams.external.CopyTrack(tr);
	beams.external.SetLocation(false);
	scatteredBeams.push_back(beams.external);

	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);

	if (IsTerminalAct(beams.internal))
	{
		if (!beams.internal.isInside)
		{
			ReleaseBeam(beams.internal);
		}
	}
	else
	{
		PushBeamToTree(beams.internal);
	}

#ifdef _CHECK_ENERGY_BALANCE
	ComputeFacetEnergy(facet->in_normal, beams.external);
#endif
}

void ScatteringConvex::SplitOriginalBeam(std::vector<Beam> &externalBeams)
{
	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, 0, m_particle->nElems, m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		m_splitting.SetBeams(*facet);
		ComputeOpticalBeamParams(facet, m_originalBeam);
		PushBeamsToBuffer(facet, m_splitting.beams, externalBeams);
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	FindVisibleFacets(beam, m_lightChecker, 0, m_particle->nElems, facets);
}
