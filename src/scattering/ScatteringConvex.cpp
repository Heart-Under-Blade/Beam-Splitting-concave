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

		Beam inBeam, outBeam;
		SplitLightToBeams(facet, inBeam, outBeam);

		auto newId = RecomputeTrackId(0, facet->index);

		outBeam.id = newId;
		outBeam.facet = facet;
		outBeam.nActs = 0;
		outBeams.push_back(outBeam);

		inBeam.id = newId;
		PushBeamToTree(inBeam, facet, 0, Location::In);

#ifdef _CHECK_ENERGY_BALANCE
		ComputeFacetEnergy(facet->in_normal, outBeam);
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

			Beam inBeam;
			bool isIncident = SplitSecondaryBeams(beam, facet, inBeam, outBeams);

			if (!isIncident)
			{
				continue;
			}

			inBeam.id = RecomputeTrackId(beam.id, i);
			inBeam.locations = beam.locations;
			PushBeamToTree(inBeam, facet, beam.nActs+1, Location::In);
		}
	}
}

bool ScatteringConvex::SplitSecondaryBeams(Beam &incidentBeam, Facet *facet,
										   Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.direction;

	// ext. normal uses in this calculating
	const Point3f &normal = facet->ex_normal;
	m_splitting.ComputeCosA(normal, incidentDir);

	if (!m_splitting.IsIncident())
	{
		return false;
	}

	Intersect(facet, incidentBeam, outBeam);

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
			outBeam.id = RecomputeTrackId(incidentBeam.id, facet->index);
			outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
			outBeam.facet = facet;
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
		outBeam.id = RecomputeTrackId(incidentBeam.id, facet->index);
		double path = m_splitting.ComputeOutgoingOpticalPath(outBeam); // добираем оптический путь
		outBeam.opticalPath += path;
		outBeam.facet = facet;
		outBeams.push_back(outBeam);
	}

	return true;
}
