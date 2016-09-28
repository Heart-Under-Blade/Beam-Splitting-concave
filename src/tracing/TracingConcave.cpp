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
		double cosIncident = DotProduct(m_startBeam.direction, extNormal);

		Beam inBeam, outBeam;
		SetBeamsParamsExternal(facetIndex, cosIncident, inBeam, outBeam);

		if (i == 0) // facet in front of beam
		{
			SetBeamShapesByFacet(facetIndex, inBeam, outBeam);
		}
		else // probably shadowed facet
		{
			SetBeamShapesByClipping(facetIndices, i, true, inBeam, outBeam);
		}

		if (m_isOpticalPath)
		{
			Point3f center = inBeam.Center();

			inBeam.D = DotProduct(-inBeam.direction, center);
			inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_startBeam.direction, center);

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

void TracingConcave::SelectVisibleFacets(const BeamInfoConcave &info, int *facetIndices,
										 int &indicesNumber)
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

void TracingConcave::SetPolygonByFacet(const Point3f *facet, int size, Paths &polygon)
{
	for (int i = 0; i < size; ++i)
	{
		const Point3f &p = facet[i];

		polygon[0] << IntPoint((int)std::round(p.cx * MULTI_INDEX),
							   (int)std::round(p.cy * MULTI_INDEX),
							   (int)std::round(p.cz * MULTI_INDEX));
	}
}

void TracingConcave::CutShadowsOutOfFacet(int *facetIndices, int facetCount, const Point3f &normal,
										  Paths &result)
{
	Paths origin(1);
	int facetIndex = facetIndices[facetCount]; // current facet index
	Point3f *facet = m_particle->facets[facetIndex];
	int size = m_particle->vertexNums[facetIndex];
	SetPolygonByFacet(facet, size, origin);

	for (int i = 0; i < facetCount; ++i)
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

void TracingConcave::SetBeamShapeByPolygon(Beam &beam, const Paths &result)
{
	int vertexNum = result.size();
	beam.shapeSize = vertexNum;

	for (const IntPoint &p : result[0])
	{
		Point3f tmp(p.X / MULTI_INDEX,
					p.Y / MULTI_INDEX,
					p.Z / MULTI_INDEX);

		beam.shape[--vertexNum] = tmp;
	}
}

void TracingConcave::SetBeamShapesByClipping(int *facetIndices, int facetCount, bool isExternal,
											 Beam &inBeam, Beam &outBeam)
{
	Point3f &normal = (isExternal) ? m_particle->externalNormals[facetCount]
									 : m_particle->normals[facetCount];

	Paths result;
	CutShadowsOutOfFacet(facetIndices, facetCount, normal, result);

	SetBeamShapeByPolygon(inBeam, result);
	SetBeamShapeByPolygon(outBeam, result);
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
	Point3f *facet = m_particle->facets[facetIndex];

	for (int i = 0; i < m_particle->vertexNums[facetIndex]; ++i)
	{
		// measure dist
		Point3f point;
		ProjectPointToFacet(facet[i], -beamDir, beamDir, point);
		double newDist = sqrt(Norm(point - facet[i]));

		if (newDist < dist) // choose minimum with previews
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

void TracingConcave::CutBeamShape(const Beam &outBeam, Beam &incidentBeam)
{
	Paths origin(1);
	SetPolygonByFacet(incidentBeam.shape, incidentBeam.shapeSize, origin);

	Paths cutter, result;
	SetPolygonByFacet(outBeam.shape, outBeam.shapeSize, cutter);

	Clipper clipper;
	clipper.AddPaths(origin, ptSubject, true);
	clipper.AddPaths(cutter, ptClip, true);
	clipper.Execute(ctXor /* excepted difference of polygons */,
					result, pftNonZero, pftNonZero);

	SetBeamShapeByPolygon(incidentBeam, result);
}

void TracingConcave::TraceInternalReflections(BeamInfoConcave *tree, int size,
											  std::vector<Beam> &outBeams)
{
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

		incidentDir.D_PARAM = m_particle->normals[info.facetIndex].D_PARAM;
		SelectVisibleFacets(info, facetIndices, indexCount);

		SortFacets(indexCount, incidentDir, facetIndices);

		for (int i = 0; i < m_particle->facetNum; ++i)
		{
			int facetIndex = facetIndices[i];

			if (facetIndex == info.facetIndex)
			{
				continue;
			}

			Beam inBeam;

			Beam outBeam;
			const Point3f &normal = m_particle->externalNormals[facetIndex];
			double cosIncident = DotProduct(normal, incidentDir);

			bool isOk = Intersect(facetIndex, incidentBeam, outBeam);

			if (!isOk)
			{
				if (info.isExternal) // beam has left the particle
				{
					outBeams.push_back(outBeam);
				}

				continue;
			}

			if (info.isExternal)
			{
				CutBeamShape(outBeam, incidentBeam);
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

			tree[size++] = BeamInfoConcave(inBeam, facetIndex, info.dept+1, false);
		}

		if (info.isExternal) // посылам обрезанный всеми гранями внешний пучок на сферу
		{
			outBeams.push_back(incidentBeam);
		}
	}
}


void TracingConcave::findDividePoint(const Point3f *polygon, int size,
									 int i0, int i1, const Point3f &normal,
									 Point3f &x, int &nextPointIndex) const
{
	__m128 res;

	__m128 _a0 = _mm_load_ps(polygon[i0]);
	__m128 _a1 = _mm_load_ps(polygon[i1]);

	bool isOk = false;
	int j0 = size-1;

	for (int i = 0; i < size; ++i)
	{
		__m128 _b0 = _mm_load_ps(polygon[j0]);
		__m128 _b1 = _mm_load_ps(polygon[i]);

		res = computeIntersection_i(_a0, _a1, _b0, _b1, normal, isOk); // TODO: написать опт. вариант (с уже готовыми векторами v1, v2)

		if (isOk)
		{
			x.point[0] = res[0];
			x.point[1] = res[1];
			x.point[2] = res[2];
			nextPointIndex = i;
			return;
		}
	}

	assert(false && "Divide point is not found");
}

void TracingConcave::FillSubpolygon(int begin, int end,
									const Point3f *polygon, int size,
									std::vector<Point3f> &subpolygon) const
{
	for (int j = begin; j != end; ++j)
	{
		if (j == size)
		{
			j = 0;
		}

		subpolygon.push_back(polygon[j]);
	}
}

void TracingConcave::DivideConcavePolygon(const Point3f *polygon, int size,
										  const Point3f &normal,
										  Polygons &polygons) const
{
	std::vector<Point3f> subPolygon;

	int baseI = size-1;
	Point3f v1 = polygon[0] - polygon[baseI];

	subPolygon.push_back(polygon[baseI]);

	for (int i = 1; i < size; ++i)
	{
		subPolygon.push_back(polygon[i]);

		Point3f v2 = polygon[i] - polygon[baseI];
		Point3f res;
		CrossProduct(v2, v1, res);
		double dir = DotProduct(res, normal);

		if (dir < -EPS_INTERSECTION) // cavity is here
		{
			Point3f x;
			int next;

			findDividePoint(polygon, size, baseI, i-1, normal, x, next);

			std::vector<Point3f> pol1;
			pol1.push_back(polygon[i-1]);
			pol1.push_back(x);
			FillSubpolygon(next, i-1, polygon, size, pol1);

			std::vector<Point3f> pol2;
			pol2.push_back(polygon[next]);
			pol2.push_back(x);
			FillSubpolygon(i-1, next, polygon, size, pol2);



			DivideConcavePolygon();
		}

		++baseI;
		v1 = v2; // TODO: опт. с индексами
	}
}

double TracingConcave::BeamCrossSection(const Beam &beam) const
{
	const double Eps = 1e7*DBL_EPSILON;

	Point3f normal;
	Point3f p1 = beam.shape[1] - beam.shape[0];
	Point3f p2 = beam.shape[2] - beam.shape[0];
	CrossProduct(p2, p1, normal);

	normal.cx = fabs(normal.cx);
	normal.cy = fabs(normal.cy);
	normal.cz = fabs(normal.cz);

	double e = fabs(DotProduct(normal, beam.direction));

	if (e < Eps)
	{
		return 0;
	}

	DivideConcavePolygon(beam.shape, beam.shapeSize, normal);

	double square = 0;
	{
		const Point3f &basePoint = beam.shape[0];
		Point3f v1 = beam.shape[1] - basePoint;

		for (int i = 2; i < beam.shapeSize; ++i)
		{
			Point3f v2 = beam.shape[i] - basePoint;
			Point3f res;
			CrossProduct(v2, v1, res);
			square += sqrt(Norm(res));
			v1 = v2;
		}

		if (square < 0)
		{	/// TODO: для опт. узнать в какую сторону ориентированы точки в пучке
			square *= (-1);
		}

		square /= 2.0;
	}

	double n = sqrt(Norm(normal));
	return (e*square) / n; // TODO: опт.
}
