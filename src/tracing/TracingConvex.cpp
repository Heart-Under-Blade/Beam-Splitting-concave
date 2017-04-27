#include "TracingConvex.h"

TracingConvex::TracingConvex(Particle *particle, const Point3f &incidentBeamDir, bool isOpticalPath,
							 const Point3f &polarizationBasis, int interReflectionNumber)
	: Tracing(particle, incidentBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
}

void TracingConvex::SplitBeamByParticle(double beta, double gamma, std::vector<Beam> &outBeams)
{
	m_particle->Rotate(beta, gamma, 0);

	m_lightSurfaceArea = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int facetID = 0; facetID < m_particle->m_facetNum; ++facetID)
	{
		const Point3f &extNormal = m_particle->facets[facetID].ex_normal;
		double cosIN = DotProduct(m_waveFront.direction, extNormal);

		if (cosIN >= EPS_M_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;

//		if (facetID == 3) // DEB
//			int ggg = 0;
		TraceFirstBeam(facetID, inBeam, outBeam);

		outBeam.lastFacetID = facetID;
		outBeam.level = 0;
		SetBeamID(outBeam);
		outBeams.push_back(outBeam);
		PushBeamToTree(inBeam, facetID, 0);

		if (m_isArea)
		{
			CalcLigthSurfaceArea(facetID, outBeam);
		}
	}

	TraceInternalReflections(outBeams);
}

void TracingConvex::SplitBeamByParticle(double, double, const std::vector<std::vector<int>> &/*tracks*/, std::vector<Beam> &)
{
}

void TracingConvex::TraceInternalReflections(std::vector<Beam> &outBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_beamTree[--m_treeSize];

		if (IsTerminalBeam(beam))
		{
			continue;
		}

		for (int facetID = 0; facetID < m_particle->m_facetNum; ++facetID)
		{
			if (facetID == beam.lastFacetID)
			{
				continue;
			}

			Beam inBeam;

			try
			{
				TraceSecondaryBeams(beam, facetID, inBeam, outBeams);
			}
			catch (const std::exception &)
			{
				continue;
			}

			inBeam.id = beam.id;
			PushBeamToTree(inBeam, facetID, beam.level+1);
		}
	}
}
