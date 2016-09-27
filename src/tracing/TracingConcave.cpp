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
	m_startBeam.shape[0] = m_startBeam.direction * m_particle->halfHeight*2;
	++m_startBeam.shapeSize;
	m_startBeam.direction.D_PARAM = DotProduct(m_startBeam.shape[0], -m_startBeam.direction);
}

void TracingConcave::SplitBeamByParticle(std::vector<Beam> &outBeams,
										 double &lightSurfaceSquare)
{
	BeamInfoConcave tree[MAX_BEAM_DEPT]; /// beam info tree (based on stack)
	int size = 0;

	/// first extermal beam

	int facetIndices[m_particle->facetNum];
	int indexCount = 0;

	BeamInfoConcave info;
	info.beam = m_startBeam;
	info.isExternal = true;

	SelectVisibleFacets(info, facetIndices, indexCount);

	SortFacets(indexCount, m_startBeam.direction, facetIndices);

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
		tree[size++] = BeamInfoConcave(inBeam, facetIndex, 0, false);
		tree[size++] = BeamInfoConcave(outBeam, facetIndex, 0, true);
		lightSurfaceSquare += outBeam.Square() * cosIncident;
	}

	TraceInternalReflections(tree, size, outBeams);
}

void TracingConcave::SelectVisibleFacets(const BeamInfoConcave &info, int *facetIndices, int &indicesNumber)
{
	Point3f *normals = (info.isExternal) ? m_particle->externalNormals
										 : m_particle->normals;

	for (int i = 0; i < m_particle->facetNum; ++i)
	{
		double cosIncident = DotProduct(info.beam.direction, normals[i]);

		if (cosIncident >= EPS_COS89) /// beam incidents to this facet
		{
			Point3f vectorToFacet = m_particle->facets[i][0] - info.beam.shape[0];
			double cosFacets = DotProduct(info.beam.direction, vectorToFacet);

			if (cosFacets < EPS_COS89) /// facet is behind the begin of beam
			{
				facetIndices[indicesNumber] = i;
				++indicesNumber;
			}
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

void TracingConcave::CutShadowsOutOfFacet(int *facetIndices, int currIndex, Paths &result)
{
	Paths origin(1);
	int facetIndex = facetIndices[currIndex];
	SetPolygonByFacet(facetIndex, origin);

	Point3f &normal = (m_isBeamInParticle) ? m_particle->normals[facetIndex]
											: m_particle->externalNormals[facetIndex];

	for (int i = 0; i < currIndex; ++i)
	{
		Paths cutter(1);
		ProjectFacetToFacet(facetIndices[i], normal, cutter);

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
	CutShadowsOutOfFacet(facetIndices, previewFacetCount, result);

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

void TracingConcave::ProjectPointToFacet(const Point3f &point, const Point3f &direction,
										 const Point3f &facetNormal, Point3f &projection)
{
	double t = DotProduct(point, facetNormal);
	t = t + facetNormal.D_PARAM;
	double dp = DotProduct(direction, facetNormal);
	t = t/dp;
	projection = point - (direction * t);
}

void TracingConcave::ProjectFacetToFacet(int a_index, const Point3f &normal,
										 Paths &projection)
{
	Point3f *a_facet = m_particle->facets[a_index];
	int &size = m_particle->vertexNums[a_index];

	for (int i = 0; i < size; ++i)
	{
		Point3f p;
		ProjectPointToFacet(a_facet[i], -normal, normal, p);

		projection[0] << IntPoint((int)std::round(p.cx * MULTI_INDEX),
								  (int)std::round(p.cy * MULTI_INDEX),
								  (int)std::round(p.cz * MULTI_INDEX));
	}
}

double TracingConcave::MeasureMinDistanceToFacet(int facetIndex, const Point3f &beamDir)
{
	double dist = FLT_MAX;
	const Point3f &facet = m_particle->facets[facetIndex];

	for (int i = 0; i < m_particle->vertexNums[facetIndex]; ++i)
	{
		// measure dist
		Point3f point;
		ProjectPointToFacet(facet[i], -beamDir, beamDir, point);
		double newDist = sqrt(Norm(point - facet[i]));

		if (newDist < dist) // choose minimum with preview
		{
			dist = newDist;
		}
	}

	return dist;
}

void TracingConcave::SortFacets(int number, const Point3f &beamDir, int *facetIndices)
{
	float distances[number];

	for (int i = 0; i < number; ++i)
	{
		distances[i] = MeasureMinDistanceToFacet(facetIndices[i], beamDir);
	}

	int left = 0;
	int rigth = number - 1;

	int stack[MAX_VERTEX_NUM*2];
	int size = 0;

	stack[size++] = left;
	stack[size++] = rigth;

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
			stack[size++] = i;
			stack[size++] = rigth;
		}

		if (left < j)
		{
			stack[size++] = left;
			stack[size++] = j;
		}

		if (size == 0)
		{
			break;
		}

		rigth = stack[--size];
		left = stack[--size];
	}
}

void TracingConcave::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
										 std::vector<Beam> &outBeams)
{
}

void TracingConcave::TraceInternalReflections(BeamInfo *tree, int size,
											  std::vector<Beam> &outBeams)
{
	m_isBeamInParticle = true;

	while (size != 0)
	{
		BeamInfoConcave info = tree[--size];

		if (isEnough(info))
		{
			continue;
		}

		Beam &incidentBeam = info.beam;
		Point3f &incidentDir = incidentBeam.direction;

		int facetIndices[m_particle->facetNum];
		int indexCount = 0;

		incidentDir.D_PARAM = m_particle->facets[info.facetIndex];
		SelectVisibleFacets(incidentBeam, facetIndices, indexCount);

		SortFacets(indexCount, incidentDir, facetIndices);

		for (int i = 0; i < m_particle->facetNum; ++i)
		{
			int facetIndex = facetIndices[i];

			if (facetIndex == info.facetIndex)
			{
				continue;
			}

			Beam inBeam;

			try
			{
				Beam outBeam;
				const Point3f &normal = m_particle->externalNormals[facetIndex];
				double cosIncident = DotProduct(normal, incidentDir);

				bool isOk = Intersect(facetIndex, incidentBeam, outBeam);

				if (!isOk)
				{
					throw std::exception();
				}

				inBeam = outBeam;
				double cosInc_sqr = cosIncident * cosIncident;

				double Nr;
				{
					double &re = m_particle->refrI_sqr_re;
					double &im = m_particle->refrI_sqr_im;
					Nr = (re + sqrt(re*re + im/cosInc_sqr))/2.0;
				}

				const complex &refrIndex = m_particle->refractionIndex;

				if (cosIncident >= EPS_COS0) /// case of the normal incidence
				{
					complex temp;

					temp = (2.0*refrIndex)/(1.0 + refrIndex); /// NOTE: для опт. выделить эту константу
					SetBeam(outBeam, incidentBeam, incidentDir, incidentBeam.e,
							temp, temp);

					temp = (1.0 - refrIndex)/(1.0 + refrIndex);
					SetBeam(inBeam, incidentBeam, -incidentDir, incidentBeam.e,
							temp, -temp);

					if (m_isOpticalPath)
					{
						CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
					}

					outBeams.push_back(outBeam);
				}
				else
				{
					Point3f r0 = incidentDir/cosIncident - normal;

					Point3f reflDir = r0 - normal;
					Normalize(reflDir);
					inBeam.direction = reflDir;

					Point3f scatteringNormal;
					CrossProduct(normal, incidentDir, scatteringNormal);
					Normalize(scatteringNormal);
					inBeam.e = scatteringNormal;

					incidentBeam.RotatePlane(scatteringNormal);

					double s = 1.0/(Nr*cosInc_sqr) - Norm(r0);
					complex tmp0 = refrIndex*cosIncident;

					if (s > DBL_EPSILON)
					{
						Point3f refrDir = r0/sqrt(s) + normal;
						Normalize(refrDir);

						double cosRefr = DotProduct(normal, refrDir);

						complex tmp1 = refrIndex*cosRefr;
						complex tmp = 2.0*tmp0;
						complex Tv0 = tmp1 + cosIncident;
						complex Th0 = tmp0 + cosRefr;

						SetBeam(outBeam, incidentBeam, refrDir, scatteringNormal,
								tmp/Tv0, tmp/Th0);

						complex Tv = (cosIncident - tmp1)/Tv0;
						complex Th = (tmp0 - cosRefr)/Th0;
						SetBeam(inBeam, incidentBeam, reflDir, scatteringNormal, Tv, Th);

						if (m_isOpticalPath)
						{
							CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
						}

						outBeams.push_back(outBeam);
					}
					else /// case of the complete internal reflection
					{
						const double bf = Nr*(1.0 - cosInc_sqr) - 1.0;

						double im = (bf > 0) ? sqrt(bf)
											 : 0;

						const complex sq(0, im);
						complex tmp = refrIndex*sq;
						complex Rv = (cosIncident - tmp)/(tmp + cosIncident);
						complex Rh = (tmp0 - sq)/(tmp0 + sq);

						SetBeam(inBeam, incidentBeam, reflDir, scatteringNormal, Rv, Rh);

						if (m_isOpticalPath)
						{
							Point3f center = outBeam.Center();
							inBeam.D = DotProduct(-center, inBeam.direction);

							double temp = DotProduct(incidentDir, center);

							inBeam.opticalPath = incidentBeam.opticalPath
									+ sqrt(Nr)*fabs(temp + incidentBeam.D);
						}
					}
				}
			}
			catch (const std::exception &)
			{
				continue;
			}

			tree[size++] = BeamInfo{inBeam, facetIndex, info.dept+1};
		}
	}
}
