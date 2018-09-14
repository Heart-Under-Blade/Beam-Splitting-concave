#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::ScatterLight(std::vector<Beam> &outBeams)
{
	m_incidentEnergy = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int i = 0; i < m_particle->nElems; ++i)
	{
		Facet *facet = m_particle->GetActualFacet(i);
		m_splitting.ComputeCosA(m_incidentDir, facet->in_normal);

		if (!m_splitting.IsIncident()) /// beam is not incident to this facet
		{
			continue;
		}

		SplitLightToBeams(facet);

		auto newId = RecomputeTrackId(0, facet->index);

		m_splitting.outBeam.id = newId;
		m_splitting.outBeam.facet = facet;
		m_splitting.outBeam.nActs = 0;
		m_splitting.outBeam.location = Location::Out;
		outBeams.push_back(m_splitting.outBeam);

		m_splitting.inBeam.id = newId;
		PushBeamToTree(m_splitting.inBeam, facet, 0, Location::In);

#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facet->in_normal, m_splitting.outBeam);
#endif
	}

	TraceInternalBeams(outBeams);
}

void ScatteringConvex::ScatterLight(const std::vector<std::vector<int>> &/*tracks*/,
									std::vector<Beam> &)
{
}

void ScatteringConvex::TraceInternalBeams(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_propagatingBeams[--m_treeSize];

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

			SplitSecondaryBeams(beam, facet, outBeams);
		}
	}
}

void ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, Facet *facet,
										   std::vector<Beam> &outBeams)
{
	// ext. normal uses in this calculating
	m_splitting.ComputeCosA(facet->ex_normal, incidentBeam.direction);

	if (!m_splitting.IsIncident())
	{
		return;
	}

	Polygon beamShape;
	bool isIntersected = IncidentBeamToFacet(facet, incidentBeam, beamShape);

	if (isIntersected)
	{
		m_splitting.SetBeams(beamShape);

		bool hasOutBeam = SetOpticalBeamParams(facet, incidentBeam);
		auto newId = RecomputeTrackId(incidentBeam.id, facet->index);

		if (hasOutBeam)
		{
			m_splitting.outBeam.nActs = incidentBeam.nActs + 1;
			m_splitting.outBeam.id = newId;
			m_splitting.outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(m_splitting.outBeam); // добираем оптический путь
			m_splitting.outBeam.facet = facet;
			m_splitting.outBeam.location = Location::Out;
			outBeams.push_back(m_splitting.outBeam);
		}

		PushBeamToTree(m_splitting.inBeam, incidentBeam, newId, facet, Location::In);
	}
}
