#include "Tracing.h"

#include "Mueller.hpp" /// TODO: delete
#include <float.h>
#include <assert.h>

#include "vector_lib.h"
#include <iostream>

#define EPS_COS89	1.7453292519943295769148298069306e-10	//cos(89.99999999)
#define EPS_COS0	0.99999999998254670756866631966593		//1- cos(89.99999999)

/// TODO: алгоритм для заранее заданных траекторий

Tracing::Tracing(Particle *p_particle, bool p_isOpticalPath, const Point3f &p_incidentDir,
				 const Point3f &p_polarizationBasis, int p_interReflectionNumber)
{
	particle = p_particle;
	isOpticalPath = p_isOpticalPath;
	incidentDirection = p_incidentDir;
	polarizationBasis = p_polarizationBasis;
	interReflectionNumber = p_interReflectionNumber;
}

double Tracing::TraceParticle(double beta, double gamma)
{
	particle->Rotate(beta, gamma, 0); /// TODO: для опт. заменить на ф-цию с двумя арг

	BeamInfo tree[MAX_BEAM_DEPT]; /// beam info tree (based on stack)
	int dept = 0;

	double square = 0;

	const complex refrIndex = particle->refractionIndex;

	/// first extermal beam
	for (int facetIndex = 0; facetIndex < particle->facetNum; ++facetIndex)
	{
		const Point3f &normal = particle->normals[facetIndex];
		const Point3f &extNormal = particle->externalNormals[facetIndex];

		double cosIncident = DotProduct(incidentDirection, extNormal);

		if (cosIncident < EPS_COS89) /// beam is not incident to this facet
		{
			continue;
		}

		Beam inBeam, outBeam;

		if (cosIncident >= EPS_COS0) /// case of the normal incidence
		{
			inBeam.JMatrix.m11 = 2.0/(refrIndex + 1.0);
			inBeam.JMatrix.m22 = inBeam.JMatrix.m11;
			inBeam.e = polarizationBasis;
			inBeam.direction = -incidentDirection;

			outBeam.JMatrix.m11 = (refrIndex - 1.0)/(refrIndex + 1.0);
			outBeam.JMatrix.m22 = -outBeam.JMatrix.m11;
			outBeam.e = polarizationBasis;
			outBeam.direction = incidentDirection;
		}
		else
		{
			/// rotate polarization plane
			{
				Point3f newBasis;
				CrossProduct(normal, incidentDirection, newBasis);
				Normalize(newBasis);
				inBeam.e = polarizationBasis;
				inBeam.direction = -incidentDirection;
				inBeam.RotatePlane(newBasis);
			}

			Point3f refrDir, reflDir;
			SplitIncidentDirection(incidentDirection, cosIncident, facetIndex,
								   reflDir, refrDir);

			double cosRefr = DotProduct(normal, reflDir);

			complex Tv00 = refrIndex*cosIncident;
			complex Th00 = refrIndex*cosRefr;

			complex Tv0 = Tv00 + cosRefr;
			complex Th0 = Th00 + cosIncident;

			complex Tv = (Tv00 - cosRefr)/Tv0;
			complex Th = (cosIncident - Th00)/Th0;
			SetBeam(outBeam, inBeam, refrDir, inBeam.e, Tv, Th);


//			using namespace std;
//			cout << endl << "! "
//				 << real(outBeam.JMatrix.m11) << ", "
//				 << real(outBeam.JMatrix.m12) << ", "
//				 << real(outBeam.JMatrix.m21) << ", "
//				 << real(outBeam.JMatrix.m22) << endl;
//			cout << endl << "! "
//				 << real(inBeam.JMatrix.m11) << ", "
//				 << real(inBeam.JMatrix.m12) << ", "
//				 << real(inBeam.JMatrix.m21) << ", "
//				 << real(inBeam.JMatrix.m22) << endl;

			double cosInc2 = (2.0*cosIncident);
			SetBeam(inBeam, inBeam, reflDir, inBeam.e,
					cosInc2/Tv0, cosInc2/Th0);
		}

		SetBeamShapesByFacet(facetIndex, inBeam, outBeam);

		if (isOpticalPath)
		{
			Point3f center = inBeam.Center();

			inBeam.D = DotProduct(-inBeam.direction, center);
			inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(incidentDirection, center);

			outBeam.D = DotProduct(-outBeam.direction, center);
			outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
		}

		tree[dept++] = BeamInfo{inBeam, facetIndex, 0};

		square += outBeam.Square()*cosIncident;

//		outBeam.RotateSpherical(incidentDirection, polarizationBasis);

		outcomingBeams.push_back(OutBeam(outBeam, track, trackSize));
	}

	TraceInternal(tree, dept);
/// TODO: проверить track'и
	trackSize = 0;
	return square;
}

void Tracing::ClearTracing()
{
	outcomingBeams.clear();
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

	return (j_norm < LOW_ENERGY_LEVEL) || (info.dept >= interReflectionNumber);
}

inline void Tracing::changeTrack(int &lastBeamDept, const BeamInfo &info)
{
	if (info.dept < lastBeamDept)
	{
		trackSize -= (lastBeamDept - info.dept);
	}
	else if (info.dept == lastBeamDept)
	{
		--trackSize;
	}

	track[trackSize++] = info.lastFacetIndex;
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

void Tracing::TraceInternal(BeamInfo tree[], int treeDept)
{
	++trackSize;
	const complex &refrIndex = particle->refractionIndex;

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
		Point3f &incidentDir = incidentBeam.direction;
		lastBeamDept = info.dept;

		for (int facetIndex = 0; facetIndex < particle->facetNum; ++facetIndex)
		{
			if (facetIndex == lastFacetIndex)
			{
				continue;
			}

			const Point3f &normal = particle->externalNormals[facetIndex];
			double cosIncident = DotProduct(normal, incidentDir);

			if (cosIncident < EPS_COS89)
			{
				continue;
			}

			Beam outBeam;
			bool isOk = Intersect(facetIndex, incidentBeam, outBeam);

			if (!isOk)
			{
				continue;
			}

			Beam inBeam = outBeam;

			double cosInc_sqr = cosIncident*cosIncident;
			double re = particle->refrI_sqr_re;
			double im = particle->refrI_sqr_im;
			double Nr = (re + sqrt(re*re + im/cosInc_sqr))/2.0;

			if (cosIncident >= EPS_COS0) /// case of the normal incidence
			{
				complex temp;

				temp = 2.0*refrIndex/(1.0 + refrIndex); /// NOTE: для опт. выделить эту константу
				SetBeam(outBeam, incidentBeam, incidentDir, incidentBeam.e,
						temp, temp);

				temp = (1.0 - refrIndex)/(1.0 + refrIndex);
				SetBeam(inBeam, incidentBeam, -incidentDir, incidentBeam.e,
						temp, -temp);

				if (isOpticalPath)
				{
					CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
				}

				/// TODO: обеспечить единую ориентацию вершин у пучков
				/// или переориентировать их тут
//				InvertBeamPointOrder(outBeam, inBeam);

				outcomingBeams.push_back(OutBeam(outBeam, track, trackSize));
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

					if (isOpticalPath)
					{
						CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
					}

					/// TODO: обеспечить единую ориентацию вершин у пучков
					/// или переориентировать их тут
//					InvertBeamPointOrder(outBeam, inBeam);

					outcomingBeams.push_back(OutBeam(outBeam, track, trackSize));
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

					if (isOpticalPath)
					{
						Point3f center = outBeam.Center();
						inBeam.D = DotProduct(-center, inBeam.direction);

						double temp = DotProduct(incidentDir, center);

						inBeam.opticalPath = incidentBeam.opticalPath
								+ sqrt(Nr)*fabs(temp + incidentBeam.D);
					}
				}
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
	__m128 _direction = _mm_setr_ps(dir.X, dir.Y, dir.Z, 0.0);

	__m128 _d_param = _mm_set_ps1(particle->normals[facetIndex].D_PARAM);
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
		__m128 _point = _mm_setr_ps(p.X, p.Y, p.Z, 0.0);
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
	__m128 _output_points[MAX_VERTEX_NUM];
	__m128 *_output_ptr = _output_points;

	int outputSize = originBeam.shapeSize;

	const Point3f &normal = particle->externalNormals[facetIndex];
	__m128 _normal_to_facet = _mm_setr_ps(normal.X, normal.Y, normal.Z, 0.0);

	const Point3f &normal2 = particle->normals[facetIndex];
	__m128 _normal_to_facet2 = _mm_setr_ps(normal2.X, normal2.Y, normal2.Z, 0.0);

	bool isProjected = ProjectToFacetPlane(originBeam, _output_points,
										   _normal_to_facet2, facetIndex);

	if (!isProjected)
	{
		return false;
	}

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	int facetSize = particle->vertexNums[facetIndex];

	__m128 _p1, _p2; /// vertices of facet
	__m128 _s_point, _e_point; /// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = particle->facets[facetIndex][facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = particle->facets[facetIndex][i];
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
			p.X = _output_points[i][0];
			p.Y = _output_points[i][1];
			p.Z = _output_points[i][2];
			outputBeam.AddVertex(p);
		}

		p0 = _output_points[i];
	}
}

void Tracing::SplitIncidentDirection(const Point3f &incidentDir, double cosIncident, int normalIndex,
									 Point3f &reflectionDir, Point3f &refractionDir) const
{
	const Point3f &normal = particle->externalNormals[normalIndex];
	const double &re = particle->refrI_sqr_re;
	const double &im = particle->refrI_sqr_im;

	Point3f tmp0 = normal - incidentDir/cosIncident;

	refractionDir = normal + tmp0;
	Normalize(refractionDir);

	double cosI_sqr = cosIncident*cosIncident;
	double tmp1 = re + cosI_sqr - 1.0;

	tmp1 = (im - FLT_EPSILON/*TODO*/ < 0) ? tmp1
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

	int vertexNum = particle->vertexNums[facetIndex];
	inBeam.shapeSize = vertexNum;
	outBeam.shapeSize = vertexNum;
	--vertexNum;

	for (int i = 0; i <= vertexNum; ++i)
	{
		inBeam.shape[i] = particle->facets[facetIndex][vertexNum-i];
		outBeam.shape[i] = particle->facets[facetIndex][vertexNum-i];
	}
}
