#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::SplitLightToBeams(std::vector<Beam> &scatteredBeams)
{
	Array<Facet*> facets;
	FindVisibleFacets(m_originBeam, m_lightChecker, 0, m_particle->nElems, facets);

	for (int i = 0; i < facets.nElems; ++i)
	{
		Facet *facet = facets.elems[i];

		m_splitting.SetBeams(*facet);
		ComputeOpticalBeamParams(facet, m_originBeam);

		auto newId = RecomputeTrackId(0, facet->index);
		m_splitting.outBeam.id = newId;
		m_splitting.outBeam.SetTracingParams(facet, 0, false);
		scatteredBeams.push_back(m_splitting.outBeam);

		m_splitting.inBeam.id = newId;
		PushBeamToTree(m_splitting.inBeam, facet, 0, true);

#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facet->in_normal, m_splitting.outBeam);
#endif
	}
}

void ScatteringConvex::SplitBeams(std::vector<Beam> &scatteredBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_tracingBeams[--m_treeSize];

		if (IsTerminalAct(beam))
		{
			continue;
		}
		else
		{
			SplitBeamByFacets(beam, scatteredBeams);
		}
	}
}

void ScatteringConvex::SplitBeamByFacets(Beam &beam, std::vector<Beam> &scatteredBeams)
{
	Array<Facet*> facets;
	FindVisibleFacets(beam, m_lightChecker, 0, m_particle->nElems, facets);

	for (int i = 0; i < facets.nElems; ++i)
	{
		Facet *facet = facets.elems[i];

		Polygon beamShape;
		bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, beam.isInside,
														   beam.direction, beamShape);
		if (isIntersected)
		{
			// ext. normal uses in this calculating
			m_splitting.SetBeams(beamShape);
			bool hasOutBeam = ComputeOpticalBeamParams(facet, beam);

			auto newId = RecomputeTrackId(beam.id, facet->index);

			if (hasOutBeam)
			{
				m_splitting.outBeam.SetTracingParams(facet, beam.act+1, false);
				m_splitting.outBeam.id = newId;
				m_splitting.outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(m_splitting.outBeam); // добираем оптический путь
				scatteredBeams.push_back(m_splitting.outBeam);
			}

			PushBeamToTree(m_splitting.inBeam, beam, newId, facet, true);
		}
	}
}
