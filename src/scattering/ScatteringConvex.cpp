#include "ScatteringConvex.h"

ScatteringConvex::ScatteringConvex(Particle *particle, Light *incidentLight,
								   bool isOpticalPath, int nActs)
	: Scattering(particle, incidentLight, isOpticalPath, nActs)
{
}

void ScatteringConvex::ScatterLight(double beta, double gamma, std::vector<Beam> &outBeams)
{
	m_particle->Rotate(beta, gamma, 0);

	m_incommingEnergy = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int facetID = 0; facetID < m_particle->nFacets; ++facetID)
	{
		const Point3f &extNormal = m_facets[facetID].ex_normal;
		double cosIN = DotProduct(m_incidentDir, extNormal);

		if (cosIN >= EPS_M_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		SplitLightToBeams(facetID, inBeam, outBeam);

		outBeam.lastFacetId = facetID;
		outBeam.level = 0;
		ComputeBeamId(outBeam);
		outBeams.push_back(outBeam);
		PushBeamToTree(inBeam, facetID, 0);

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

			inBeam.trackId = beam.trackId;
			PushBeamToTree(inBeam, id, beam.level+1);
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
	double cosA = DotProduct(normal, incidentDir);

	if (cosA < EPS_COS_90) /// beam is not incident to this facet
	{
		return false;
	}

	bool isOk = Intersect(facetID, incidentBeam, outBeam);

	if (!isOk)
	{
		return false;
	}

	inBeam = outBeam;

	if (cosA < EPS_COS_00)
	{	// regular incidence
		bool isTrivialIncidence;
		SetRegularIncidenceBeamParams(cosA, normal, incidentBeam,
									  inBeam, outBeam, isTrivialIncidence);
		if (isTrivialIncidence)
		{
			outBeam.trackId = incidentBeam.trackId;
			outBeam.lastFacetId = facetID;
			outBeam.level = incidentBeam.level + 1;
			ComputeBeamId(outBeam);
			outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.frontPosition); // добираем оптический путь
			outBeams.push_back(outBeam);
		}
	}
	else
	{	// normal incidence
		SetNormalIncidenceBeamParams(cosA, incidentBeam, inBeam, outBeam);

		outBeam.trackId = incidentBeam.trackId;
		outBeam.lastFacetId = facetID;
		outBeam.level = incidentBeam.level + 1;
		ComputeBeamId(outBeam);
		outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.frontPosition); // добираем оптический путь
		outBeams.push_back(outBeam);
	}

	return true;
}
