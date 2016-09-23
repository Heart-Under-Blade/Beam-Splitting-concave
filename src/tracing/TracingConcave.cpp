#include "TracingConcave.h"

#include <assert.h>
#include <float.h>
#include <tgmath.h>

#define MULTI_INDEX 10000000 // index for poligon's clip operations

using namespace ClipperLib;

TracingConcave::TracingConcave(Particle *particle, const Point3f &startBeamDir,
							   bool isOpticalPath, const Point3f &polarizationBasis,
							   int interReflectionNumber)
	: Tracing(particle, startBeamDir, isOpticalPath, polarizationBasis, interReflectionNumber)
{
	Point3f facetPoint = m_startBeamDirection * m_particle->halfHeight*2;
	m_startBeamDParam = DotProduct(facetPoint, -m_startBeamDirection);
}

void TracingConcave::SplitBeamByParticle(std::vector<OutBeam> &outBeams,
										 double &lightSurfaceSquare)
{
	BeamInfo tree[MAX_BEAM_DEPT]; /// beam info tree (based on stack)
	int treeSize = 0;

	/// first extermal beam

	int facetIndices[m_particle->facetNum];
	int indexCount = 0;
	SelectVisibleFacets(facetIndices, indexCount);

	SortFacets(facetIndices, indexCount);

	for (int i = 0; i < indexCount; ++i)
	{
		int facetIndex = facetIndices[i];
		const Point3f &extNormal = m_particle->externalNormals[facetIndex];
		double cosIncident = DotProduct(m_startBeamDirection, extNormal);

		Beam inBeam, outBeam;
		SetBeamsParamsExternal(facetIndex, cosIncident, inBeam, outBeam);

		if (i == 0) // facet in front of beam
		{
			SetBeamShapesByFacet(facetIndex, inBeam, outBeam);
		}
		else // probably shadowed facet
		{
			SetBeamShapesByClipping(facetIndices, i, inBeam, outBeam);
		}

		if (m_isOpticalPath)
		{
			Point3f center = inBeam.Center();

			inBeam.D = DotProduct(-inBeam.direction, center);
			inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_startBeamDirection, center);

			outBeam.D = DotProduct(-outBeam.direction, center);
			outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
		}

//		outBeams.push_back(OutBeam(outBeam, m_track, m_trackSize));
		tree[treeSize++] = BeamInfo{inBeam, facetIndex, 0};
		tree[treeSize++] = BeamInfo{outBeam, facetIndex, 0};
		lightSurfaceSquare += outBeam.Square()*cosIncident;
	}
}

void TracingConcave::SelectVisibleFacets(int *facetIndices, int &indicesNumber)
{
	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		const Point3f &extNormal = m_particle->externalNormals[i];
		double cosIncident = DotProduct(m_startBeamDirection, extNormal);

		if (cosIncident >= EPS_COS89) /// beam is not incident to this facet
		{
			facetIndices[indicesNumber] = i;
			++indicesNumber;
		}
	}
}

void TracingConcave::SetPolygonByFacet(int facetIndex, Paths &polygon)
{
	for (int i = 0; i < m_particle->vertexNums[facetIndex]; ++i)
	{
		Point3f &p = m_particle->facets[facetIndex][i];

		polygon[0] << IntPoint((int)std::round(p.cx * MULTI_INDEX),
							   (int)std::round(p.cy * MULTI_INDEX),
							   (int)std::round(p.cz * MULTI_INDEX));
	}
}

void TracingConcave::CutShadesOutOfFacet(int *facetIndices, int previewFacetCount, Paths &result)
{
	Paths origin(1);
	SetPolygonByFacet(facetIndices[previewFacetCount], origin);

	for (int i = 0; i < previewFacetCount; ++i)
	{
		Paths cutter(1);
		SetPolygonByFacet(facetIndices[i], cutter);

		Clipper clipper;
		clipper.AddPaths(origin, ptSubject, true);
		clipper.AddPaths(cutter, ptClip, true);
		clipper.Execute(ctXor /* excepted difference of polygons */,
						result, pftNonZero, pftNonZero);

		if (result.size() < MIN_VERTEX_NUM)
		{
			continue; // clipping operation is not succeed
		}

		origin = result;
	}

	assert(result.size() < MIN_VERTEX_NUM); // polygon must be closed
}

void TracingConcave::SetBeamShapesByClipping(int *facetIndices, int previewFacetCount,
											 Beam &inBeam, Beam &outBeam)
{
	Paths result;
	CutShadesOutOfFacet(facetIndices, previewFacetCount, result);

	// fill beam shapes
	int vertexNum = result.size();
	inBeam.shapeSize = vertexNum;
	outBeam.shapeSize = vertexNum;

	for (const IntPoint &p : result[0])
	{
		Point3f tmp(p.X / MULTI_INDEX,
					p.Y / MULTI_INDEX,
					p.Z / MULTI_INDEX);

		inBeam.shape[--vertexNum] = tmp;
		outBeam.shape[--vertexNum] = tmp;
	}
}

double TracingConcave::MeasureMinDistanceToFacet(int facetIndex)
{
	double dist = FLT_MAX;

	for (int i = 0; i < m_particle->vertexNums[facetIndex]; ++i)
	{
		// measure dist
		Point3f &p = m_particle->facets[facetIndex][i];
		double t = DotProduct(p, -m_startBeamDirection);
		t = t + m_startBeamDParam;
		double dp = Norm(m_startBeamDirection);
		t = t/dp;
		Point3f point = m_startBeamDirection * t;
		point = p - point;
		double newDist = sqrt(Norm(point - p));

		if (newDist < dist) // choose min
		{
			dist = newDist;
		}
	}

	return dist;
}

void TracingConcave::SortFacets(int *facetIndices, int number)
{
	float distances[number];

	for (int i = 0; i < number; ++i)
	{
		distances[i] = MeasureMinDistanceToFacet(facetIndices[i]);
	}

	int left = 0;
	int rigth = number - 1;

	int stack[MAX_VERTEX_NUM*2];
	int stackSize = 0;

	stack[stackSize++] = left;
	stack[stackSize++] = rigth;

	while (true)
	{
		int base = distances[(left + rigth)/2];

		int i = left;
		int j = rigth;

		while (i <= j)
		{
			while (distances[i] < base)
			{
				++i;
			}

			while (distances[j] > base)
			{
				--j;
			}

			if (i <= j)	// exchange elems
			{
				float temp_d = distances[i];
				distances[i] = distances[j];
				distances[j] = temp_d;

				int temp_v = facetIndices[i];
				facetIndices[i] = facetIndices[j];
				facetIndices[j] = temp_v;

				++i;
				--j;
			}
		}

		if (i < rigth)
		{
			stack[stackSize++] = i;
			stack[stackSize++] = rigth;
		}

		if (left < j)
		{
			stack[stackSize++] = left;
			stack[stackSize++] = j;
		}

		if (stackSize == 0)
		{
			break;
		}

		rigth = stack[--stackSize];
		left = stack[--stackSize];
	}
}
