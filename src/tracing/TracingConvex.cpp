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
	for (int facetId = 0; facetId < m_particle->facetNum; ++facetId)
	{
		const Point3f &extNormal = m_particle->facets[facetId].ex_normal;
		double cosIN = DotProduct(m_waveFront.direction, extNormal);

		if (cosIN < EPS_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		SplitExternalBeamByFacet(facetId, inBeam, outBeam);

		outBeams.push_back(outBeam);
		m_beamTree[m_treeSize] = inBeam;
		m_beamTree[m_treeSize].facetID = facetId;
		m_beamTree[m_treeSize].level = 0;
		++m_treeSize;

		if (m_isArea)
		{
			CalcLigthSurfaceArea(facetId, outBeam);
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

		if (isTerminalBeam(beam))
		{
			continue;
		}

		for (int facetIndex = 0; facetIndex < m_particle->facetNum; ++facetIndex)
		{
			if (facetIndex == beam.facetID)
			{
				continue;
			}

			Beam inBeam;

			try
			{
				SplitInternalBeamByFacet(beam, facetIndex, inBeam, outBeams);
			}
			catch (const std::exception &)
			{
				continue;
			}

			m_beamTree[m_treeSize] = inBeam;
			m_beamTree[m_treeSize].facetID = facetIndex;
			m_beamTree[m_treeSize].level = beam.level+1;
			++m_treeSize;
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
