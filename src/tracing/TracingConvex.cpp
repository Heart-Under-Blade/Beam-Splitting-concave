#include "TracingConvex.h"

TracingConvex::TracingConvex(Particle *particle, const Point3f &incidentBeamDir, bool isOpticalPath,
							 const Point3f &polarizationBasis, int interReflectionNumber)
	: Tracing(particle, incidentBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
}

void TracingConvex::SplitBeamByParticle(std::vector<Beam> &outBeams)
{
	m_lightSurfaceArea = 0;
	m_treeSize = 0;

	/// first extermal beam
	for (int facetID = 0; facetID < m_particle->facetNum; ++facetID)
	{
		const Point3f &extNormal = m_particle->facets[facetID].ex_normal;
		double cosIN = DotProduct(m_waveFront.direction, extNormal);

		if (cosIN >= EPS_M_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		TraceFirstBeam(facetID, inBeam, outBeam);

		outBeam.facetID = facetID;
		outBeam.level = 0;
		SetBeamId(outBeam);
		outBeams.push_back(outBeam);
		PushBeamToTree(inBeam, facetID, 0);

		if (m_isArea)
		{
			CalcLigthSurfaceArea(facetID, outBeam);
		}
	}

	TraceInternalReflections(outBeams);
}

void TracingConvex::SplitBeamByParticle(const std::vector<std::vector<int>> &/*tracks*/,
										std::vector<Beam> &/*outBeams*/)
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

		for (int facetID = 0; facetID < m_particle->facetNum; ++facetID)
		{
			if (facetID == beam.facetID)
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


double TracingConvex::BeamCrossSection(const Beam &beam) const
{
	const double Eps = 1e7*DBL_EPSILON;

	Point3f normal;
	Point3f p1 = beam.polygon.arr[1] - beam.polygon.arr[0];
	Point3f p2 = beam.polygon.arr[2] - beam.polygon.arr[0];
	CrossProduct(p1, p2, normal);

	double e = fabs(DotProduct(normal, beam.direction));

	if (e < Eps)
	{
		return 0;
	}

	double square = 0;
	{
		const Point3f &basePoint = beam.polygon.arr[0];
		Point3f p1 = beam.polygon.arr[1] - basePoint;

		for (int i = 2; i < beam.polygon.size; ++i)
		{
			Point3f p2 = beam.polygon.arr[i] - basePoint;
			Point3f res;
			CrossProduct(p1, p2, res);
			square += sqrt(Norm(res));
			p1 = p2;
		}

		square /= 2.0;
	}

	double n = sqrt(Norm(normal));
	return (e*square) / n;
}
