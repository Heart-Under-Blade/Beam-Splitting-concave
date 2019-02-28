#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::ScatterLight(double beta, double gamma,
									std::vector<Beam> &outBeams)
{
	m_particle->Rotate(beta, gamma, 0);

	m_incidentEnergy = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int facetID = 0; facetID < m_particle->nFacets; ++facetID)
	{
		const Point3f &inNormal = m_facets[facetID].in_normal;
		m_splitting.ComputeCosA(m_incidentDir, inNormal);

		if (!m_splitting.IsIncident()) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		SplitLightToBeams(facetID, inBeam, outBeam);

		auto newId = RecomputeTrackId(0, facetID);

		outBeam.id = newId;
		outBeam.lastFacetId = facetID;
		outBeam.nActs = 0;
		outBeams.push_back(outBeam);

		inBeam.id = newId;
		PushBeamToTree(inBeam, facetID, 0, Location::In);

#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facetID, outBeam);
#endif
	}

	TraceInternalBeams(outBeams);
}

void ScatteringConvex::ScatterLight(double, double, const std::vector<std::vector<int>> &/*tracks*/, std::vector<Beam> &)
{
}

void ScatteringConvex::TraceInternalBeams(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];

		if (IsTerminalAct(beam))
		{
			continue;
		}

		for (int id = 0; id < m_particle->nFacets; ++id)
		{
			if (id == beam.lastFacetId)
			{
				continue;
			}

			Beam inBeam;
			bool isIncident = SplitSecondaryBeams(beam, id, inBeam, outBeams);

			if (!isIncident)
			{
				continue;
			}

			inBeam.id = RecomputeTrackId(beam.id, id);
			inBeam.locations = beam.locations;
			PushBeamToTree(inBeam, id, beam.nActs+1, Location::In);
		}
	}
}

bool ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, int facetID,
										   Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.direction;

	// ext. normal uses in this calculating
	const Point3f &normal = m_facets[facetID].ex_normal;
	m_splitting.ComputeCosA(normal, incidentDir);

	if (!m_splitting.IsIncident())
	{
		return false;
	}

	Intersect(facetID, incidentBeam, outBeam);

	if (outBeam.nVertices < MIN_VERTEX_NUM)
	{
		return false;
	}

	inBeam = outBeam;
	m_splitting.ComputeSplittingParams(incidentBeam.direction, normal);

	if (!m_splitting.IsNormalIncidence())
	{	// regular incidence
		ComputePolarisationParams(incidentBeam.direction, normal, incidentBeam);

		if (!m_splitting.IsCompleteReflection())
		{
			m_splitting.ComputeRegularBeamsParams(normal, incidentBeam,
												  inBeam, outBeam);
			outBeam.nActs = incidentBeam.nActs + 1;
			outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);
			outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
			outBeam.lastFacetId = facetID;
			outBeams.push_back(outBeam);
		}
		else // complete internal reflection incidence
		{
			m_splitting.ComputeCRBeamParams(normal, incidentBeam, inBeam);
		}
	}
	else
	{	// normal incidence
		m_splitting.ComputeNormalBeamParams(incidentBeam, inBeam, outBeam);

		outBeam.nActs = incidentBeam.nActs + 1;
		outBeam.id = RecomputeTrackId(incidentBeam.id, facetID);
		double path = m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
		outBeam.opticalPath += path;
		outBeam.lastFacetId = facetID;
		outBeams.push_back(outBeam);
	}

	return true;
}
