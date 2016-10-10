#include "TracingConvex.h"

TracingConvex::TracingConvex(Particle *particle, const Point3f &incidentBeamDir, bool isOpticalPath,
							 const Point3f &polarizationBasis, int interReflectionNumber)
	: Tracing(particle, incidentBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
}

void TracingConvex::SplitBeamByParticle(std::vector<Beam> &outBeams,
										double &lightSurfaceSquare)
{
	/// TODO: отделить функцию высчитывания площади осв. поверхности
	lightSurfaceSquare = 0;

	BeamInfo tree[MAX_BEAM_DEPT]; /// beam info tree (based on stack)
	int treeSize = 0;

	/// first extermal beam
	for (int facetIndex = 0; facetIndex < m_particle->facetNum; ++facetIndex)
	{
		const Point3f &extNormal = m_particle->externalNormals[facetIndex];
		double cosIncident = DotProduct(m_startBeam.direction, extNormal);

		if (cosIncident < EPS_COS_89) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;
		SplitExternalBeamByFacet(facetIndex, cosIncident, inBeam, outBeam);

		outBeams.push_back(outBeam);
		tree[treeSize++] = BeamInfo{inBeam, facetIndex, 0};
//		lightSurfaceSquare += outBeam.Square()*cosIncident;
	}

	TraceInternalReflections(tree, treeSize, outBeams);
}

void TracingConvex::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
										std::vector<Beam> &outBeams)
{

}

void TracingConvex::TraceInternalReflections(BeamInfo *tree, int treesize,
											 std::vector<Beam> &outBeams)
{
	while (treesize != 0)
	{
		BeamInfo info = tree[--treesize];

		if (isEnough(info))
		{
			continue;
		}

		for (int facetIndex = 0; facetIndex < m_particle->facetNum; ++facetIndex)
		{
			if (facetIndex == info.facetId)
			{
				continue;
			}

			Beam inBeam;

			try
			{
				SplitInternalBeamByFacet(info.beam, facetIndex, inBeam, outBeams);
			}
			catch (const std::exception &)
			{
				continue;
			}

			tree[treesize++] = BeamInfo{inBeam, facetIndex, info.dept+1};
		}
	}
}


double TracingConvex::BeamCrossSection(const Beam &beam) const
{
	const double Eps = 1e7*DBL_EPSILON;

	Point3f normal;
	Point3f p1 = beam.shape[1] - beam.shape[0];
	Point3f p2 = beam.shape[2] - beam.shape[0];
	CrossProduct(p1, p2, normal);

	double e = fabs(DotProduct(normal, beam.direction));

	if (e < Eps)
	{
		return 0;
	}

	double square = 0;
	{
		const Point3f &basePoint = beam.shape[0];
		Point3f p1 = beam.shape[1] - basePoint;

		for (int i = 2; i < beam.shapeSize; ++i)
		{
			Point3f p2 = beam.shape[i] - basePoint;
			Point3f res;
			CrossProduct(p1, p2, res);
			square += sqrt(Norm(res));
			p1 = p2;
		}

		if (square < 0)
		{	/// OPT: узнать в какую сторону ориентированы точки в пучке
			square *= (-1);
		}

		square /= 2.0;
	}

	double n = sqrt(Norm(normal));
	return (e*square) / n; // OPT:
}
