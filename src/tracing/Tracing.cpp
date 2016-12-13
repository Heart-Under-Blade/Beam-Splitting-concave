#include "Tracing.h"

//#include "Mueller.hpp"
#include <float.h>
#include <assert.h>

#include "macro.h"
#include "vector_lib.h"
//#include <iostream>

#define NORM_CEIL	FLT_EPSILON + 1

Tracing::Tracing(Particle *particle, const Point3f &startBeamDir, bool isOpticalPath,
				 const Point3f &polarizationBasis, int interReflectionNumber)
{
	LOG_ASSERT(startBeamDir.cx <= NORM_CEIL
		   && startBeamDir.cy <= NORM_CEIL
		   && startBeamDir.cz <= NORM_CEIL
		   && "Direction of the start beam is not normalized.");

	m_particle = particle;
	m_isOpticalPath = isOpticalPath;
	m_polarizationBasis = polarizationBasis;
	m_interReflectionNumber = interReflectionNumber;
	m_startBeam.direction = startBeamDir;
}

void Tracing::RotateParticle(double beta, double gamma)
{
	m_particle->Rotate(beta, gamma, 0);
}

void Tracing::SetBeamsParamsExternal(int facetIndex, double cosIncident,
									 Beam &inBeam, Beam &outBeam)
{
	const complex &refrIndex = m_particle->refractionIndex;
	const Point3f &facetNormal = m_particle->normals[facetIndex];
	const Point3f &startDir = m_startBeam.direction;

	if (cosIncident < EPS_COS_00)
	{
		/// rotate polarization plane
		{
			Point3f newBasis;
			CrossProduct(facetNormal, -startDir, newBasis);
			Normalize(newBasis);
			inBeam.e = m_polarizationBasis;
			inBeam.direction = startDir;
			inBeam.RotatePlane(newBasis);
		}

		Point3f refrDir, reflDir;
		SplitBeamDirection(startDir, cosIncident, -facetNormal, reflDir, refrDir);

		double cosRelr = DotProduct(facetNormal, reflDir);

		complex Tv00 = refrIndex*cosIncident;
		complex Th00 = refrIndex*cosRelr;

		complex Tv0 = Tv00 + cosRelr;
		complex Th0 = Th00 + cosIncident;

		complex Tv = (Tv00 - cosRelr)/Tv0;
		complex Th = (cosIncident - Th00)/Th0;
		SetBeam(outBeam, inBeam, refrDir, inBeam.e, Tv, Th);

		double cosInc2 = (2.0*cosIncident);
		SetBeam(inBeam, inBeam, reflDir, inBeam.e,
				cosInc2/Tv0, cosInc2/Th0);
	}
	else /// normal incidence
	{
		inBeam.JMatrix.m11 = 2.0/(refrIndex + 1.0);
		inBeam.JMatrix.m22 = inBeam.JMatrix.m11;
		inBeam.e = m_polarizationBasis;
		inBeam.direction = -startDir;

		outBeam.JMatrix.m11 = (refrIndex - 1.0)/(refrIndex + 1.0);
		outBeam.JMatrix.m22 = -outBeam.JMatrix.m11;
		outBeam.e = m_polarizationBasis;
		outBeam.direction = startDir;
	}
}

void Tracing::SplitExternalBeamByFacet(int facetIndex, double cosIncident,
									   Beam &inBeam, Beam &outBeam)
{
	SetBeamsParamsExternal(facetIndex, cosIncident, inBeam, outBeam);

	SetBeamByFacet(facetIndex, inBeam);
	SetBeamByFacet(facetIndex, outBeam);

	if (m_isOpticalPath)
	{
		Point3f center = CenterOfPolygon(inBeam.polygon, inBeam.size);

		inBeam.D = DotProduct(-inBeam.direction, center);
		inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_startBeam.direction, center);

		outBeam.D = DotProduct(-outBeam.direction, center);
		outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
	}
}

// TODO: пофиксить
void Tracing::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
								  std::vector<Beam> &outBeams)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		int facetIndex = tracks.at(i).at(0);
		const Point3f &extNormal = m_particle->externalNormals[facetIndex];

		double cosIncident = DotProduct(m_startBeam.direction, extNormal);

		if (cosIncident < EPS_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		std::vector<Beam> outBuff;
		Beam incidentBeam;

		/// first incident beam
		{
			Beam outBeam;

			SplitExternalBeamByFacet(facetIndex, cosIncident, incidentBeam, outBeam);

			outBuff.push_back(outBeam);
		}

		unsigned int size = tracks.at(i).size();

		try /// internal beams
		{
			for (unsigned int j = 1; j < size; ++j)
			{
				facetIndex = tracks.at(i).at(j);

				Beam inBeam;
				SplitInternalBeamByFacet(incidentBeam, facetIndex, inBeam, outBuff);

				incidentBeam = inBeam;
			}
		}
		catch (const std::exception &)
		{
			continue;
		}

		outBeams.push_back(outBuff.back());
	}
}

void Tracing::CalcOpticalPathInternal(double Nr, const Beam &incidentBeam,
									  Beam &outBeam, Beam &inBeam) const
{
	double coef = (incidentBeam.isExternal) ? 1 : sqrt(Nr);
	Point3f center = CenterOfPolygon(outBeam.polygon, outBeam.size);

	outBeam.D = DotProduct(-outBeam.direction, center);

	double temp = DotProduct(incidentBeam.direction, center);
	outBeam.opticalPath = incidentBeam.opticalPath
			+ coef*fabs(temp + incidentBeam.D);

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = outBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + inBeam.D);
}

bool Tracing::isEnough(const Beam &beam)
{
	double j_norm = beam.JMatrix.Norm();
	return (j_norm < LOW_ENERGY_LEVEL) || (beam.level >= m_interReflectionNumber);
}

void Tracing::InvertBeamShapeOrder(Beam &outBeam, const Beam &inBeam)
{
	Point3f center = CenterOfPolygon(outBeam.polygon, outBeam.size);
	Point3f p0 = outBeam.polygon[0] - center;
	Point3f p1 = outBeam.polygon[1] - center;
	Point3f normal;
	CrossProduct(p0, p1, normal);
	double cosAngle = DotProduct(normal, outBeam.direction);

	if (cosAngle < 0.0)
	{
		for (int i = 0; i < inBeam.size; ++i)
		{
			outBeam.polygon[i] = inBeam.polygon[(inBeam.size-1)-i];
		}
	}
}

void Tracing::SplitInternalBeamByFacet(Beam &incidentBeam, int facetIndex,
									   Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;

	const Point3f &incidentDir = incidentBeam.direction;
	const complex &refrIndex = m_particle->refractionIndex;

	const Point3f &normal = m_particle->externalNormals[facetIndex]; /// ext. normal uses in this calculating
	double cosIncident = DotProduct(normal, incidentDir);

	if (cosIncident < EPS_COS_90) /// beam is not incident to this facet
	{
		throw std::exception();
	}

	bool isOk = Intersect(facetIndex, incidentBeam, outBeam);

	if (!isOk)
	{
		throw std::exception();
	}

	inBeam = outBeam;
	double cosInc_sqr = cosIncident*cosIncident;

	double Nr;
	{
		double re = m_particle->refrI_coef_re;
		double im = m_particle->refrI_coef_im;
		Nr = (re + sqrt(re*re + im/cosInc_sqr))/2.0;
	}

	if (cosIncident >= EPS_COS_00) /// case of the normal incidence
	{
		complex temp;

		temp = (2.0*refrIndex)/(1.0 + refrIndex);
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
				Point3f center = CenterOfPolygon(outBeam.polygon, outBeam.size);
				inBeam.D = DotProduct(-center, inBeam.direction);

				double temp = DotProduct(incidentDir, center);

				inBeam.opticalPath = incidentBeam.opticalPath
						+ sqrt(Nr)*fabs(temp + incidentBeam.D);
			}
		}
	}
}

void Tracing::SetBeam(Beam &beam, const Beam &other, const Point3f &dir, const Point3f &e,
					  const complex &coef1, const complex &coef2) const
{
	beam.MulJMatrix(other, coef1, coef2);

	beam.direction = dir;
	beam.e = e;
}

/// Projection of beam to facet
bool Tracing::ProjectToFacetPlane(const Point3f *polygon, int size, const Point3f &dir,
								  const Point3f &normal, __m128 *_projection) const
{
	__m128 _normal = _mm_setr_ps(normal.cx, normal.cy, normal.cz, 0.0);
	__m128 _direction = _mm_setr_ps(dir.cx, dir.cy, dir.cz, 0.0);

	__m128 _d_param = _mm_set_ps1(normal.d_param);
	__m128 _dp0 = _mm_dp_ps(_direction, _normal, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_PROJECTION)
	{
		return false; /// beam is parallel to facet
	}

	for (int i = 0; i < size; ++i)
	{
		const Point3f &p = polygon[i];
		__m128 _point = _mm_setr_ps(p.cx, p.cy, p.cz, 0.0);
		__m128 _dp1 = _mm_dp_ps(_point, _normal, MASK_FULL);
		__m128 _add = _mm_add_ps(_dp1, _d_param);
		__m128 _t = _mm_div_ps(_add, _dp0);
		__m128 _mul = _mm_mul_ps(_t, _direction);

		_projection[i] = _mm_sub_ps(_point, _mul);
	}

	return true;
}

/// NOTE: вершины пучка и грани должны быть ориентированы в одном направлении
bool Tracing::Intersect(int facetId, const Beam &beam, Beam &intersection) const
{
	__m128 _output_points[MAX_VERTEX_NUM];

	const Point3f &normal = m_particle->normals[facetId];
	bool isProjected = ProjectToFacetPlane(beam.polygon, beam.size, beam.direction,
										   normal, _output_points);
	if (!isProjected)
	{
		return false;
	}

	__m128 _normal_to_facet = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
	__m128 *_output_ptr = _output_points;
	int outputSize = beam.size;

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	int facetSize = m_particle->vertexNums[facetId];

	__m128 _p1, _p2; /// vertices of facet
	__m128 _s_point, _e_point; /// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = m_particle->facets[facetId][facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = m_particle->facets[facetId][i];
		_p2 = _mm_load_ps(p2.point);

		bufferSize = outputSize;
		outputSize = 0;

		__m128 *_temp = _output_ptr;
		_output_ptr = _buffer_ptr;
		_buffer_ptr = _temp;

		_s_point = _buffer_ptr[bufferSize-1];
		isInsideS = is_inside_i(_s_point, _p1, _p2, _normal_to_facet);

		bool isIntersected;

		for (int j = 0; j < bufferSize; ++j)
		{
			_e_point = _buffer_ptr[j];
			isInsideE = is_inside_i(_e_point, _p1, _p2, _normal_to_facet);

			if (isInsideE)
			{
				if (!isInsideS)
				{
					__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2, _normal_to_facet, isIntersected);

					if (isIntersected)
					{
						_output_ptr[outputSize++] = x;
					}
				}

				_output_ptr[outputSize++] = _e_point;
			}
			else if (isInsideS)
			{
				__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2, _normal_to_facet, isIntersected);

				if (isIntersected)
				{
					_output_ptr[outputSize++] = x;
				}
			}

			_s_point = _e_point;
			isInsideS = isInsideE;
		}
	}

	SetOutputBeam(_output_ptr, outputSize, intersection);
	return intersection.size >= MIN_VERTEX_NUM;
}

void Tracing::SetOutputBeam(__m128 *_output_points, int outputSize, Beam &outputBeam) const
{
	Point3f p;

	__m128 eps = _mm_load_ps1(&EPS_INTERSECTION);
	__m128 sign_mask = _mm_set1_ps(-0.f);

	__m128 p0 = _output_points[outputSize-1];

	for (int i = 0; i < outputSize; ++i)
	{
		__m128 difference = _mm_sub_ps(_output_points[i], p0);
		__m128 abs = _mm_andnot_ps(sign_mask, difference);
		__m128 cmp = _mm_cmplt_ps(eps, abs);

		int res = _mm_movemask_ps(cmp) & 0b111;

		if (res != 0)
		{
			p.cx = _output_points[i][0];
			p.cy = _output_points[i][1];
			p.cz = _output_points[i][2];
			outputBeam.AddVertex(p);
		}

		p0 = _output_points[i];
	}
}

void Tracing::SplitBeamDirection(const Point3f &incidentDir, double cosIncident,
								 const Point3f &normal,
								 Point3f &refleDir, Point3f &refraDir) const
{
	const double &re = m_particle->refrI_coef_re;
	const double &im = m_particle->refrI_coef_im;

	Point3f tmp0 = normal + incidentDir/cosIncident;

	refraDir = normal + tmp0;
	Normalize(refraDir);

	double cosI_sqr = cosIncident * cosIncident;
	double tmp1 = re + cosI_sqr - 1.0;

	tmp1 = (im - FLT_EPSILON < 0) ? tmp1
								  : sqrt(tmp1*tmp1 + im);

	tmp1 = (re + 1.0 - cosI_sqr + tmp1)/2.0;
	tmp1 = (tmp1/cosI_sqr) - Norm(tmp0);

	LOG_ASSERT(tmp1 > 0);
	tmp1 = sqrt(tmp1);

	refleDir = (tmp0/tmp1) - normal;
	Normalize(refleDir);
}

/** NOTE: Result beams are ordered in inverse direction */
void Tracing::SetBeamByFacet(int facetId, Beam &beam) const
{
	int vertexNum = m_particle->vertexNums[facetId];
	beam.size = vertexNum;
	--vertexNum;

	for (int i = 0; i <= vertexNum; ++i)
	{
		beam.polygon[i] = m_particle->facets[facetId][vertexNum-i];
	}
}

double Tracing::Square(const Beam &beam)
{
	double square = 0;
	const Point3f &basePoint = beam.polygon[0];
	Point3f p1 = beam.polygon[1] - basePoint;

	for (int i = 2; i < beam.size; ++i)
	{
		Point3f p2 = beam.polygon[i] - basePoint;
		Point3f res;
		CrossProduct(p1, p2, res);
		square += sqrt(Norm(res));
		p1 = p2;
	}

	return square / 2.0;
}
