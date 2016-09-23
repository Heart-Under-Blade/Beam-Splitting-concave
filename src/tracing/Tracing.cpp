#include "Tracing.h"

//#include "Mueller.hpp"
#include <float.h>
#include <assert.h>

#include "vector_lib.h"
//#include <iostream>

#define NORMALIZE_CEIL 1 + FLT_EPSILON

Tracing::Tracing(Particle *particle, const Point3f &startBeamDir, bool isOpticalPath,
				 const Point3f &polarizationBasis, int interReflectionNumber)
{
	assert(startBeamDir.cx <= NORMALIZE_CEIL
		   && startBeamDir.cx <= NORMALIZE_CEIL
		   && startBeamDir.cz <= NORMALIZE_CEIL
		   && "Direction of the start beam is not normalized.");

	m_trackSize = 0;
	m_particle = particle;
	m_isOpticalPath = isOpticalPath;
	m_polarizationBasis = polarizationBasis;
	m_interReflectionNumber = interReflectionNumber;
	m_startBeamDirection = startBeamDir;
}

void Tracing::RotateParticle(double beta, double gamma)
{
	m_particle->Rotate(beta, gamma, 0); /// TODO: для опт. заменить на ф-цию с двумя арг
}

void Tracing::SetBeamsParamsExternal(int facetIndex, double cosIncident, Beam &outBeam, Beam &inBeam)
{
	const complex &refrIndex = m_particle->refractionIndex;
	const Point3f &normal = m_particle->normals[facetIndex];

	if (cosIncident >= EPS_COS0) /// case of the normal incidence
	{
		inBeam.JMatrix.m11 = 2.0/(refrIndex + 1.0);
		inBeam.JMatrix.m22 = inBeam.JMatrix.m11;
		inBeam.e = m_polarizationBasis;
		inBeam.direction = -m_startBeamDirection;

		outBeam.JMatrix.m11 = (refrIndex - 1.0)/(refrIndex + 1.0);
		outBeam.JMatrix.m22 = -outBeam.JMatrix.m11;
		outBeam.e = m_polarizationBasis;
		outBeam.direction = m_startBeamDirection;
	}
	else
	{
		/// rotate polarization plane
		{
			Point3f newBasis;
			CrossProduct(normal, m_startBeamDirection, newBasis);
			Normalize(newBasis);
			inBeam.e = m_polarizationBasis;
			inBeam.direction = -m_startBeamDirection;
			inBeam.RotatePlane(newBasis);
		}

		Point3f refrDir, reflDir;
		SplitIncidentDirection(m_startBeamDirection, cosIncident, normal,
							   reflDir, refrDir);

		double cosRefr = DotProduct(normal, reflDir);

		complex Tv00 = refrIndex*cosIncident;
		complex Th00 = refrIndex*cosRefr;

		complex Tv0 = Tv00 + cosRefr;
		complex Th0 = Th00 + cosIncident;

		complex Tv = (Tv00 - cosRefr)/Tv0;
		complex Th = (cosIncident - Th00)/Th0;
		SetBeam(outBeam, inBeam, refrDir, inBeam.e, Tv, Th);

		double cosInc2 = (2.0*cosIncident);
		SetBeam(inBeam, inBeam, reflDir, inBeam.e,
				cosInc2/Tv0, cosInc2/Th0);
	}
}

void Tracing::SplitExternalBeamByFacet(int facetIndex, double cosIncident,
									   Beam &inBeam, Beam &outBeam)
{
	SetBeamsParamsExternal(facetIndex, cosIncident, outBeam, inBeam);

	SetBeamShapesByFacet(facetIndex, inBeam, outBeam);

	if (m_isOpticalPath)
	{
		Point3f center = inBeam.Center();

		inBeam.D = DotProduct(-inBeam.direction, center);
		inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_startBeamDirection, center);

		outBeam.D = DotProduct(-outBeam.direction, center);
		outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
	}
}

void Tracing::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
								  std::vector<OutBeam> &outBeams)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		int facetIndex = tracks.at(i).at(0);
		const Point3f &extNormal = m_particle->externalNormals[facetIndex];

		double cosIncident = DotProduct(m_startBeamDirection, extNormal);

		if (cosIncident < EPS_COS89) /// beam is not incident to this facet
		{
			continue;
		}

		std::vector<OutBeam> buff;
		Beam incidentBeam;

		/// first incident beam
		{
			Beam outBeam;

			SplitExternalBeamByFacet(facetIndex, cosIncident, incidentBeam, outBeam);

			buff.push_back(OutBeam(outBeam, m_track, m_trackSize));
		}

		unsigned int size = tracks.at(i).size();

		try /// inteernal beams
		{
			for (unsigned int j = 1; j < size; ++j)
			{
				facetIndex = tracks.at(i).at(j);

				Beam inBeam;
				SplitInternalBeamByFacet(incidentBeam, facetIndex, inBeam, buff);

				incidentBeam = inBeam;
			}
		}
		catch (const std::exception &)
		{
			continue;
		}

		outBeams.push_back(buff.back());
	}

	m_trackSize = 0;
}

void Tracing::CalcOpticalPathInternal(double Nr, const Beam &incidentBeam,
									  Beam &outBeam, Beam &inBeam) const
{
	Point3f center = outBeam.Center();

	outBeam.D = DotProduct(-outBeam.direction, center);

	double temp = DotProduct(incidentBeam.direction, center);
	outBeam.opticalPath = incidentBeam.opticalPath
			+ sqrt(Nr)*fabs(temp + incidentBeam.D);

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = outBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + inBeam.D);
}

inline bool Tracing::isEnough(const BeamInfo &info)
{
	double j_norm = info.beam.JMatrix.Norm();

	return (j_norm < LOW_ENERGY_LEVEL) || (info.dept >= m_interReflectionNumber);
}

inline void Tracing::changeTrack(int &lastBeamDept, const BeamInfo &info)
{
	if (info.dept < lastBeamDept)
	{
		m_trackSize -= (lastBeamDept - info.dept);
	}
	else if (info.dept == lastBeamDept)
	{
		--m_trackSize;
	}

	m_track[m_trackSize++] = info.lastFacetIndex;
}

void Tracing::InvertBeamPointOrder(Beam &outBeam, const Beam &inBeam)
{
	Point3f inBeamCenter = outBeam.Center();
	Point3f p0 = outBeam.shape[0] - inBeamCenter;
	Point3f p1 = outBeam.shape[1] - inBeamCenter;
	Point3f normal1;
	CrossProduct(p1, p0, normal1);
	double cosAngle = DotProduct(normal1, outBeam.direction);

	if (cosAngle < 0.0)
	{
		for (int i = 0; i < inBeam.shapeSize; ++i)
		{
			outBeam.shape[i] = inBeam.shape[(inBeam.shapeSize-1)-i];
		}
	}
}

void Tracing::SplitInternalBeamByFacet(Beam &incidentBeam, int facetIndex,
									   Beam &inBeam, std::vector<OutBeam> &outBeams)
{
	Beam outBeam;

	const Point3f &incidentDir = incidentBeam.direction;
	const complex &refrIndex = m_particle->refractionIndex;

	const Point3f &normal = m_particle->externalNormals[facetIndex];
	double cosIncident = DotProduct(normal, incidentDir);

	if (cosIncident < EPS_COS89) /// beam is not incident to this facet
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
		double re = m_particle->refrI_sqr_re;
		double im = m_particle->refrI_sqr_im;
		Nr = (re + sqrt(re*re + im/cosInc_sqr))/2.0;
	}

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

		outBeams.push_back(OutBeam(outBeam, m_track, m_trackSize));
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

			outBeams.push_back(OutBeam(outBeam, m_track, m_trackSize));
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

void Tracing::TraceInternalReflections(BeamInfo *tree, int treeDept,
									   std::vector<OutBeam> &outBeams)
{
	++m_trackSize;
	int lastBeamDept = 0;

	while (treeDept != 0)
	{
		BeamInfo info = tree[--treeDept];

		if (isEnough(info))
		{
			continue;
		}

		changeTrack(lastBeamDept, info);

		Beam &incidentBeam = info.beam;
		int &lastFacetIndex = info.lastFacetIndex;
		lastBeamDept = info.dept;

		for (int facetIndex = 0; facetIndex < m_particle->facetNum; ++facetIndex)
		{
			if (facetIndex == lastFacetIndex)
			{
				continue;
			}

			Beam inBeam;

			try
			{
				SplitInternalBeamByFacet(incidentBeam, facetIndex, inBeam, outBeams);
			}
			catch (const std::exception &)
			{
				continue;
			}

			tree[treeDept++] = BeamInfo{inBeam, facetIndex, lastBeamDept+1};
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
bool Tracing::ProjectToFacetPlane(const Beam& inputBeam, __m128 *_output_points,
								  __m128 _normal, int facetIndex) const
{
	const Point3f &dir = inputBeam.direction;
	__m128 _direction = _mm_setr_ps(dir.cx, dir.cx, dir.cz, 0.0);

	__m128 _d_param = _mm_set_ps1(m_particle->normals[facetIndex].D_PARAM);
	__m128 _dp0 = _mm_dp_ps(_direction, _normal, MASK_FULL);

	__m128 _sign_mask = _mm_set1_ps(-0.f);
	__m128 _abs_dp = _mm_andnot_ps(_sign_mask, _dp0);

	if (_abs_dp[0] < EPS_PROJECTION)
	{
		return false; /// beam is parallel to facet
	}

	for (int i = 0; i < inputBeam.shapeSize; ++i)
	{
		const Point3f &p = inputBeam.shape[i];
		__m128 _point = _mm_setr_ps(p.cx, p.cx, p.cz, 0.0);
		__m128 _dp1 = _mm_dp_ps(_point, _normal, MASK_FULL);
		__m128 _add = _mm_add_ps(_dp1, _d_param);
		__m128 _t = _mm_div_ps(_add, _dp0);
		__m128 _mul = _mm_mul_ps(_t, _direction);

		_output_points[i] = _mm_sub_ps(_point, _mul);
	}

	return true;
}

bool Tracing::Intersect(int facetIndex, const Beam &originBeam, Beam &intersectBeam) const
{
	/// TODO: попробовать заменить на Clipper
	__m128 _output_points[MAX_VERTEX_NUM];
	__m128 *_output_ptr = _output_points;

	int outputSize = originBeam.shapeSize;

	const Point3f &normal = m_particle->externalNormals[facetIndex];
	__m128 _normal_to_facet = _mm_setr_ps(normal.cx, normal.cx, normal.cz, 0.0);

	const Point3f &normal2 = m_particle->normals[facetIndex];
	__m128 _normal_to_facet2 = _mm_setr_ps(normal2.cx, normal2.cx, normal2.cz, 0.0);

	bool isProjected = ProjectToFacetPlane(originBeam, _output_points,
										   _normal_to_facet2, facetIndex);
	if (!isProjected)
	{
		return false;
	}

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	int facetSize = m_particle->vertexNums[facetIndex];

	__m128 _p1, _p2; /// vertices of facet
	__m128 _s_point, _e_point; /// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = m_particle->facets[facetIndex][facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = m_particle->facets[facetIndex][i];
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

	SetOutputBeam(_output_ptr, outputSize, intersectBeam);
	return intersectBeam.shapeSize > 2;
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
			p.cx = _output_points[i][1];
			p.cz = _output_points[i][2];
			outputBeam.AddVertex(p);
		}

		p0 = _output_points[i];
	}
}

void Tracing::SplitIncidentDirection(const Point3f &incidentDir, double cosIncident, const Point3f &normal,
									 Point3f &reflectionDir, Point3f &refractionDir) const
{
	const double &re = m_particle->refrI_sqr_re;
	const double &im = m_particle->refrI_sqr_im;

	Point3f tmp0 = normal - incidentDir/cosIncident;

	refractionDir = normal + tmp0;
	Normalize(refractionDir);

	double cosI_sqr = cosIncident*cosIncident;
	double tmp1 = re + cosI_sqr - 1.0;

	tmp1 = (im - FLT_EPSILON < 0) ? tmp1
								  : sqrt(tmp1*tmp1 + im);

	tmp1 = (re + 1.0 - cosI_sqr + tmp1)/2.0;
	tmp1 = (tmp1/cosI_sqr) - Norm(tmp0);
	tmp1 = sqrt(tmp1);

	assert(!isnan(tmp1));
	reflectionDir = (tmp0/tmp1) - normal;
	Normalize(reflectionDir);
}

void Tracing::SetBeamShapesByFacet(int facetIndex, Beam &inBeam, Beam &outBeam) const
{
	/** NOTE: Result beams are ordered in inverse direction */

	int vertexNum = m_particle->vertexNums[facetIndex];
	inBeam.shapeSize = vertexNum;
	outBeam.shapeSize = vertexNum;
	--vertexNum;

	for (int i = 0; i <= vertexNum; ++i)
	{
		inBeam.shape[i] = m_particle->facets[facetIndex][vertexNum-i];
		outBeam.shape[i] = m_particle->facets[facetIndex][vertexNum-i];
	}
}
