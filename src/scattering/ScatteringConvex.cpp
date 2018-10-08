#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::SplitLightToBeams(std::vector<Beam> &scatteredBeams)
{
	for (int i = 0; i < m_particle->nElems; ++i)
	{
		Facet *facet = m_particle->GetActualFacet(i);
		m_splitting.ComputeCosA(m_originBeam.direction, facet->in_normal);

		if (m_splitting.IsIncident()) /// beam is not incident to this facet
		{
			m_splitting.SetBeams(*facet);
			SetOpticalBeamParams(facet, m_originBeam);

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
}

void ScatteringConvex::ScatterLight(std::vector<Beam> &scatteredBeams)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
#endif
	m_treeSize = 0;

	SplitLightToBeams(scatteredBeams);
	ScatterBeams(scatteredBeams);
}

void ScatteringConvex::ScatterLight(const std::vector<std::vector<int>> &/*tracks*/,
									std::vector<Beam> &)
{
}

void ScatteringConvex::ScatterBeams(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_tracingBeams[--m_treeSize];

		if (IsTerminalAct(beam))
		{
			continue;
		}

		for (int i = 0; i < m_particle->nElems; ++i)
		{
			Facet *facet = m_particle->GetActualFacet(i);

			if (facet->index == beam.facet->index)
			{
				continue;
			}

			SplitBeamsByFacet(beam, facet, outBeams);
		}
	}
}

void ScatteringConvex::SplitBeamsByFacet(Beam &beam, Facet *facet,
										 std::vector<Beam> &outBeams)
{
	// ext. normal uses in this calculating
	m_splitting.ComputeCosA(facet->ex_normal, beam.direction);

	if (!m_splitting.IsIncident())
	{
		return;
	}

	Polygon beamShape;
	bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, beam.isInside,
													   beam.direction, beamShape);
	if (isIntersected)
	{
		m_splitting.SetBeams(beamShape);

		bool hasOutBeam = SetOpticalBeamParams(facet, beam);
		auto newId = RecomputeTrackId(beam.id, facet->index);

		if (hasOutBeam)
		{
			m_splitting.outBeam.act = beam.act + 1;
			m_splitting.outBeam.id = newId;
			m_splitting.outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(m_splitting.outBeam); // добираем оптический путь
			m_splitting.outBeam.facet = facet;
			m_splitting.outBeam.isInside = false;
			outBeams.push_back(m_splitting.outBeam);
		}

		PushBeamToTree(m_splitting.inBeam, beam, newId, facet, true);
	}
}
