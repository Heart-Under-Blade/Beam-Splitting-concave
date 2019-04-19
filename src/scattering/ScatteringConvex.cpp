#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, const Light &incidentLight,
								   int maxActNo, const complex &refractiveIndex)
	: Scattering(particle, incidentLight, maxActNo, refractiveIndex)
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

void ScatteringConvex::ScatterLight(TrackNode *trackTree,
									std::vector<Beam> &scatteredBeams)
{
	m_hasTracks = true;
	m_trackTreeNode = trackTree;
	SetWorkFacetsFromTracks();

	SplitOriginalBeam(scatteredBeams);
	// TODO: to complete
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

	if (IsFinalAct(beams.internal))
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
	m_incidentEnergy = 0;

	m_visibleFacets.nElems = 0;
	FindVisibleFacets(m_originalBeam, m_lightChecker, m_workFacets,
					  m_visibleFacets);

	for (int i = 0; i < m_visibleFacets.nElems; ++i)
	{
		Facet *facet = m_visibleFacets.elems[i];

		ComputeOpticalBeamParams(facet, m_originalBeam, *facet);
		PushBeamsToBuffer(facet, m_splitting.beams, externalBeams);
	}
}

void ScatteringConvex::SelectVisibleFacets(const Beam &beam, Array<Facet*> &facets)
{
	if (m_hasTracks)
	{
		SetWorkFacetsFromTracks();
	}

	FindVisibleFacets(beam, m_lightChecker, m_workFacets, facets);
}
