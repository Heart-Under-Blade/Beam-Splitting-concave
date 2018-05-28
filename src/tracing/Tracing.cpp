#include "Tracing.h"

#include <float.h>
#include <assert.h>

#include "macro.h"
#include "geometry_lib.h"

#ifdef _DEBUG // DEB
#include <iostream>
#include <iomanip>
#include <limits>
#endif

#define NORM_CEIL	FLT_EPSILON + 1

using namespace std;

Tracing::Tracing(Particle *particle, Light *incidentLight, bool isOpticalPath,
				 int interReflectionNumber)
	: m_incidentLight(incidentLight),
	  m_particle(particle)
{
	m_facets = m_particle->facets;

	m_isOpticalPath = isOpticalPath;
	m_interReflectionNumber = interReflectionNumber;

	m_incidentDir = m_incidentLight->direction;
	m_incidentDir.d_param = m_incidentLight->direction.d_param;

	m_polarBasis = m_incidentLight->polarizationBasis;

	m_refrIndex = m_particle->GetRefractiveIndex();
	double re = real(m_refrIndex);
	double im = imag(m_refrIndex);
	ri_coef_re = re*re - im*im;
	ri_coef_im = 4*re*re*im;
}

void Tracing::SetSloppingBeamParams_initial(const Point3f &beamDir, double cosIN,
											int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f &facetNormal = m_particle->facets[facetId].in_normal;

	RotatePolarisationPlane(beamDir, facetNormal, inBeam);

	Point3f refrDir, reflDir;
	DivideBeamDirection(beamDir, cosIN, -facetNormal, reflDir, refrDir);

	double cosReflN = DotProduct(facetNormal, reflDir);

	complex Tv00 = m_refrIndex*cosIN;
	complex Th00 = m_refrIndex*cosReflN;

	complex Tv0 = Tv00 + cosReflN;
	complex Th0 = Th00 + cosIN;

	complex Tv = (Tv00 - cosReflN)/Tv0;
	complex Th = (cosIN - Th00)/Th0;
	SetBeam(outBeam, inBeam, refrDir, Tv, Th);

	double cosInc2 = (2.0*cosIN);
	SetBeam(inBeam, inBeam, reflDir, cosInc2/Tv0, cosInc2/Th0);
#ifdef _DEBUG // DEB
	inBeam.dirs.push_back(reflDir);
	outBeam.dirs.push_back(refrDir);
#endif
}

void Tracing::SetBeamID(Beam &beam)
{
#ifdef _TRACK_ALLOW
	beam.id += (beam.lastFacetID + 1);
	beam.id *= (m_particle->facetNum + 1);
	//	AddToTrack(beam, facetId);
#endif
}

void Tracing::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
	beam.SetTracingParams(facetId, level, location);
	PushBeamToTree(beam);
}

void Tracing::PushBeamToTree(Beam &beam, int facetId, int level)
{
	beam.lastFacetID = facetId;
	beam.level = level;
	PushBeamToTree(beam);
}

void Tracing::PushBeamToTree(Beam &beam)
{
	SetBeamID(beam);
	m_beamTree[m_treeSize++] = beam;
}

void Tracing::SetBeamOpticalParams(int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f dir = m_incidentLight->direction;
	const Point3f &normal = m_facets[facetId].in_normal;
	double cosIN = DotProduct(dir, normal);

	if (cosIN < EPS_COS_00) // slopping incidence
	{
		SetSloppingBeamParams_initial(dir, cosIN, facetId, inBeam, outBeam);
	}
	else // normal incidence
	{
		inBeam.J.m11 = 2.0/(m_refrIndex + 1.0);
		inBeam.J.m22 = inBeam.J.m11;
		inBeam.light = (*m_incidentLight);

		outBeam.J.m11 = (m_refrIndex - 1.0)/(m_refrIndex + 1.0);
		outBeam.J.m22 = -outBeam.J.m11;
		outBeam.light = Light{-dir, m_incidentLight->polarizationBasis};
#ifdef _DEBUG // DEB
		inBeam.dirs.push_back(dir);
		outBeam.dirs.push_back(-dir);
#endif
	}

	if (m_isOpticalPath)
	{
		CalcOpticalPath_initial(inBeam, outBeam);
	}
}

void Tracing::RotatePolarisationPlane(const Point3f &dir, const Point3f &facetNormal,
									  Beam &beam)
{
	Point3f newBasis;
	CrossProduct(facetNormal, -dir, newBasis);
	Normalize(newBasis);
	beam.light = Light{dir, m_incidentLight->polarizationBasis};
	beam.RotatePlane(newBasis);
}

void Tracing::CalcOpticalPath_initial(Beam &inBeam, Beam &outBeam)
{
	Point3f center = inBeam.Center();
//	Point3f &center = inBeam.arr[0];

	inBeam.D = DotProduct(-inBeam.light.direction, center);
	inBeam.opticalPath = FAR_ZONE_DISTANCE + DotProduct(m_incidentDir, center);

	outBeam.D = DotProduct(-outBeam.light.direction, center);
	outBeam.opticalPath = inBeam.opticalPath;
#ifdef _DEBUG // DEB
	inBeam.ops.push_back(inBeam.opticalPath);
	outBeam.ops.push_back(inBeam.opticalPath);
#endif
}

void Tracing::TraceFirstBeam(int facetId, Beam &inBeam, Beam &outBeam)
{
	SetPolygonByFacet(facetId, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
	SetPolygonByFacet(facetId, outBeam);
	SetBeamOpticalParams(facetId, inBeam, outBeam);
}

void Tracing::CalcFacetEnergy(int facetID, const Polygon &lightedPolygon)
{
	const Point3f &normal = m_facets[facetID].in_normal;
	double cosIN = DotProduct(m_incidentDir, normal);
	m_incommingEnergy += lightedPolygon.Area() * cosIN;
}

// TODO: пофиксить
void Tracing::SplitBeamByParticle(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
								  std::vector<Beam> &outBeams)
{
//	m_particle->Rotate(beta, gamma, 0);

//	for (unsigned int i = 0; i < tracks.size(); ++i)
//	{
//		int facetId = tracks.at(i).at(0);
//		const Point3f &extNormal = m_facets[facetId].ex_normal;

//		double cosIN = DotProduct(m_incidentDir, extNormal);

//		if (cosIN < EPS_COS_90) /// beam is not incident to this facet
//		{
//			continue;
//		}

//		std::vector<Beam> outBuff;
//		Beam incidentBeam;

//		/// first incident beam
//		{
//			Beam outBeam;
//			TraceFirstBeam(facetId, incidentBeam, outBeam);
//			outBuff.push_back(outBeam);
//		}

//		unsigned int size = tracks.at(i).size();

//		try /// internal beams
//		{
//			for (unsigned int j = 1; j < size; ++j)
//			{
//				facetId = tracks.at(i).at(j);

//				Beam inBeam;
//				TraceSecondaryBeams(incidentBeam, facetId, inBeam, outBuff);

//				incidentBeam = inBeam;
//			}
//		}
//		catch (const std::exception &)
//		{
//			continue;
//		}

//		outBeams.push_back(outBuff.back());
//	}
}

void Tracing::CalcOpticalPath(double cosIN, const Beam &incidentBeam,
							  Beam &inBeam, Beam &outBeam) const
{
	Point3f center = inBeam.Center();
//	Point3f &center = inBeam.arr[0];

	// refractive index of external environment = 1
	double OP = fabs(DotProduct(incidentBeam.light.direction, center) + incidentBeam.D);

	if (incidentBeam.location == Location::In)
	{
		OP *= 1.3116/*sqrt(CalcReRI(cosIN))*/;
//		inBeam.internalOpticalPath = outBeam.internalOpticalPath = OP;
	}

	inBeam.D = DotProduct(-inBeam.light.direction, center);
	inBeam.opticalPath = incidentBeam.opticalPath + OP;

	outBeam.D = DotProduct(-outBeam.light.direction, center);
	outBeam.opticalPath = inBeam.opticalPath;

#ifdef _DEBUG // DEB
	inBeam.ops = incidentBeam.ops;
	inBeam.ops.push_back(OP);
	outBeam.ops = incidentBeam.ops;
	outBeam.ops.push_back(OP);
#endif
}

bool Tracing::IsTerminalBeam(const Beam &beam)
{
	double j_norm = beam.J.Norm(); // OPT: move to last element of comparison
	return (j_norm < EPS_BEAM_ENERGY) || (beam.level >= m_interReflectionNumber);
}

double Tracing::CalcReRI(const double &cosIN) const
{
	double cosIN_sqr = cosIN*cosIN;
	const double &re = ri_coef_re;
	return (re + sqrt(re*re + ri_coef_im/cosIN_sqr))/2.0;
}

void Tracing::TraceSecondaryBeams(Beam &incidentBeam, int facetID,
								  Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.light.direction;

	// ext. normal uses in this calculating
	const Point3f &normal = m_particle->facets[facetID].ex_normal;
	double cosIN = DotProduct(normal, incidentDir);

	if (cosIN < EPS_COS_90) /// beam is not incident to this facet
	{
		throw std::exception();
	}

	bool isOk = Intersect(facetID, incidentBeam, outBeam);

	if (!isOk)
	{
		throw std::exception();
	}

	inBeam = outBeam;

#ifdef _DEBUG // DEB

	if (incidentBeam.lastFacetID == 0 && facetID == 5)
		int ff =0 ;
//	if (incidentBeam.id == 1278 && facetID == 2)
//		int ff =0 ;
#endif
	if (cosIN >= EPS_COS_00) /// normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);

		outBeam.id = incidentBeam.id;
		outBeam.lastFacetID = facetID;
		outBeam.level = incidentBeam.level + 1;
		SetBeamID(outBeam);
		outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.D); // добираем оптический путь
		outBeam.ops.push_back(fabs(FAR_ZONE_DISTANCE + outBeam.D));
		outBeams.push_back(outBeam);
	}
	else /// slopping incidence
	{
		bool isTrivialIncidence;
		SetSloppingIncidenceBeamParams(cosIN, normal, incidentBeam,
									   inBeam, outBeam, isTrivialIncidence);
		if (isTrivialIncidence)
		{
			outBeam.id = incidentBeam.id;
			outBeam.lastFacetID = facetID;
			outBeam.level = incidentBeam.level + 1;
			SetBeamID(outBeam);
			outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.D); // добираем оптический путь
			outBeam.ops.push_back(fabs(FAR_ZONE_DISTANCE + outBeam.D));
			outBeams.push_back(outBeam);
		}
	}
}

void Tracing::SetSloppingIncidenceBeamParams(double cosIN, const Point3f &normal,
											 Beam &incidentBeam,
											 Beam &inBeam, Beam &outBeam,
											 bool &isTrivialIncidence)
{
	const Point3f &incidentDir = incidentBeam.light.direction;
	double cosIN_sqr = cosIN*cosIN;

	Point3f scatteringNormal;
	CrossProduct(normal, incidentDir, scatteringNormal);
	Normalize(scatteringNormal);
	incidentBeam.RotatePlane(scatteringNormal);

	Point3f r0 = incidentDir/cosIN - normal;

	Point3f reflDir = r0 - normal;
	Normalize(reflDir);

	inBeam.light = Light{reflDir, scatteringNormal};
#ifdef _DEBUG // DEB
	inBeam.dirs = incidentBeam.dirs;
	inBeam.dirs.push_back(reflDir);
#endif
	double Nr = CalcReRI(cosIN);
	double s = 1.0/(Nr*cosIN_sqr) - Norm(r0);

	if (s > DBL_EPSILON) /// trivial incidence
	{
		SetTrivialIncidenceBeamParams(cosIN, Nr, normal, r0, s, incidentBeam,
									  inBeam, outBeam);
		isTrivialIncidence = true;
	}
	else /// complete internal reflection
	{
		SetCompleteReflectionBeamParams(cosIN, Nr, incidentBeam, inBeam);
		isTrivialIncidence = false;
	}
}

void Tracing::SetNormalIncidenceBeamParams(double cosIN, const Beam &incidentBeam,
										   Beam &inBeam, Beam &outBeam)
{
	const Point3f &dir = incidentBeam.light.direction;
	complex temp;

	temp = (2.0*m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(outBeam, incidentBeam, dir, temp, temp);

	temp = (1.0 - m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(inBeam, incidentBeam, -dir, temp, -temp);
#ifdef _DEBUG // DEB
	inBeam.dirs = incidentBeam.dirs;
	outBeam.dirs = incidentBeam.dirs;
	inBeam.dirs.push_back(-dir);
	outBeam.dirs.push_back(dir);
#endif
	if (m_isOpticalPath)
	{
		CalcOpticalPath(cosIN, incidentBeam, inBeam, outBeam);
	}
}

void Tracing::SetTrivialIncidenceBeamParams(double cosIN, double Nr,
											const Point3f &normal,
											Point3f r0, double s,
											const Beam &incidentBeam,
											Beam &inBeam, Beam &outBeam)
{
	Point3f refrDir = r0/sqrt(s) + normal;
	Normalize(refrDir);

	double cosRefr = DotProduct(normal, refrDir);

	complex tmp0 = m_refrIndex*cosIN;
	complex tmp1 = m_refrIndex*cosRefr;
	complex tmp = 2.0*tmp0;
	complex Tv0 = tmp1 + cosIN;
	complex Th0 = tmp0 + cosRefr;

	outBeam.SetJonesMatrix(incidentBeam, tmp/Tv0, tmp/Th0);
	outBeam.light = Light{refrDir, inBeam.light.polarizationBasis};
#ifdef _DEBUG // DEB
	outBeam.dirs = incidentBeam.dirs;
	outBeam.dirs.push_back(refrDir);
#endif
	complex Tv = (cosIN - tmp1)/Tv0;
	complex Th = (tmp0 - cosRefr)/Th0;
	inBeam.SetJonesMatrix(incidentBeam, Tv, Th);

	if (m_isOpticalPath)
	{
		CalcOpticalPath(cosIN, incidentBeam, inBeam, outBeam);
	}
}

void Tracing::SetCompleteReflectionBeamParams(double cosIN, double Nr,
											  const Beam &incidentBeam,
											  Beam &inBeam)
{
	const Point3f &incidentDir = incidentBeam.light.direction;

	double cosIN_sqr = cosIN*cosIN;
	complex tmp0 = m_refrIndex*cosIN;
	const double bf = Nr*(1.0 - cosIN_sqr) - 1.0;
	double im = (bf > 0) ? sqrt(bf) : 0;

	const complex sq(0, im);
	complex tmp = m_refrIndex*sq;
	complex Rv = (cosIN - tmp)/(tmp + cosIN);
	complex Rh = (tmp0 - sq)/(tmp0 + sq);

	inBeam.SetJonesMatrix(incidentBeam, Rv, Rh);

	if (m_isOpticalPath)
	{
		Point3f center = /*inBeam.arr[0]*/inBeam.Center();
		inBeam.D = DotProductD(-inBeam.light.direction, center);

		double temp = fabs(DotProductD(incidentDir, center) + incidentBeam.D);

		if (incidentBeam.location == Location::In)
		{
//			temp *= sqrt(Nr);
#ifdef _DEBUG // DEB
			temp *= 1.3116;
#endif
		}

		inBeam.opticalPath = incidentBeam.opticalPath + temp;
#ifdef _DEBUG // DEB
		inBeam.ops = incidentBeam.ops;
		inBeam.ops.push_back(temp);
#endif
	}
}

void Tracing::SetBeam(Beam &beam, const Beam &other, const Point3f &dir,
					  const complex &coef1, const complex &coef2) const
{
	beam.SetJonesMatrix(other, coef1, coef2);
	beam.light = Light{dir, other.light.polarizationBasis};
}

void Tracing::Difference(const Polygon &subject, const Point3f &subjNormal,
						 const Polygon &clip, const Point3f &clipNormal,
						 const Point3f &clipDir,
						 Polygon *difference, int &resultSize) const
{
	__m128 _clip[MAX_VERTEX_NUM];
	bool isProjected = ProjectToFacetPlane(clip, clipDir, subjNormal, _clip);

	if (!isProjected)
	{
		difference[resultSize++] = subject;
		return;
	}

	__m128 _clip_normal = _mm_setr_ps(clipNormal.cx, clipNormal.cy, clipNormal.cz, 0.0);

	int clipSize = clip.size;
	__m128 _diff_pol[MAX_VERTEX_NUM];

	__m128 _subject[MAX_VERTEX_NUM];
	__m128 _buffer[MAX_VERTEX_NUM];

	for (int i = 0; i < subject.size; ++i)
	{
		_subject[i] = _mm_load_ps(subject.arr[i].point);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.size;

	__m128 _first_p, _second_p;
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (int i = 0; i < clip.size; ++i)
	{
		int difSize = 0;

		__m128 _p1 = _p2;
		_p2 = _clip[i];

		__m128 *_tmp = _buff;
		_buff = _subj;
		_subj = _tmp;

		int subjSize = bufSize;
		bufSize = 0;

		_first_p = _subj[subjSize-1];
		isInFirst = is_inside_i(_first_p, _p1, _p2, _clip_normal);

		bool isIntersected;

		for (int j = 0; j < subjSize; ++j)
		{
			_second_p = _subj[j];
			isInSecond = is_inside_i(_second_p, _p1, _p2, _clip_normal);

			if (isInSecond)
			{
				if (!isInFirst)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
										   _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_diff_pol[difSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_buff[bufSize++] = _second_p;
			}
			else
			{
				if (isInFirst)
				{
					__m128 x = intersect_i(_first_p, _second_p, _p1, _p2,
										   _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_diff_pol[difSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_diff_pol[difSize++] = _second_p;
			}

			_first_p = _second_p;
			isInFirst = isInSecond;
		}

		if (difSize >= MIN_VERTEX_NUM)
		{
			Polygon resPolygon;
			SetOutputPolygon(_diff_pol, difSize, resPolygon);

			if (resPolygon.size >= MIN_VERTEX_NUM)
			{
				difference[resultSize++] = resPolygon;
			}
		}
	}
}

bool Tracing::ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
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

	for (int i = 0; i < polygon.size; ++i)
	{
		const Point3f &p = polygon.arr[i];
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
bool Tracing::Intersect(int facetID, const Beam &beam, Polygon &intersection) const
{
	__m128 _output_points[MAX_VERTEX_NUM];
	const Point3f &normal = m_facets[facetID].in_normal;

	bool isProjected = ProjectToFacetPlane(beam, beam.light.direction, normal,
										   _output_points);
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

	int facetSize = m_particle->facets[facetID].size;

	__m128 _p1, _p2; // vertices of facet
	__m128 _s_point, _e_point;	// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = m_particle->facets[facetID].arr[facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = m_particle->facets[facetID].arr[i];
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
					__m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
										   _normal_to_facet, isIntersected);
					if (isIntersected)
					{
						_output_ptr[outputSize++] = x;
					}
				}

				_output_ptr[outputSize++] = _e_point;
			}
			else if (isInsideS)
			{
				__m128 x = intersect_i(_s_point, _e_point, _p1, _p2,
									   _normal_to_facet, isIntersected);
				if (isIntersected)
				{
					_output_ptr[outputSize++] = x;
				}
			}

			_s_point = _e_point;
			isInsideS = isInsideE;
		}
	}

	SetOutputPolygon(_output_ptr, outputSize, intersection);
	return intersection.size >= MIN_VERTEX_NUM;
}

void Tracing::SetOutputPolygon(__m128 *_output_points, int outputSize,
							   Polygon &polygon) const
{
	Point3f p;

	__m128 eps = _mm_load_ps1(&EPS_MERGE);
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
			polygon.arr[polygon.size++] = p;
		}

		p0 = _output_points[i];
	}
}

void Tracing::DivideBeamDirection(const Point3f &incidentDir, double cosIN,
								  const Point3f &normal,
								  Point3f &reflDir, Point3f &refrDir) const
{
	Point3f tmp0 = normal + incidentDir/cosIN;

	refrDir = normal + tmp0;
	Normalize(refrDir);

	double cosI_sqr = cosIN * cosIN;
	double tmp1 = ri_coef_re + cosI_sqr - 1.0;

	tmp1 = (ri_coef_im < FLT_EPSILON) ? tmp1
									  : sqrt(tmp1*tmp1 + ri_coef_im);

	tmp1 = (ri_coef_re + 1.0 - cosI_sqr + tmp1)/2.0;
	tmp1 = (tmp1/cosI_sqr) - Norm(tmp0);

//	LOG_ASSERT(tmp1 > 0);
	tmp1 = sqrt(tmp1);

	reflDir = (tmp0/tmp1) - normal;
	Normalize(reflDir);
}

Point3f Tracing::ChangeBeamDirection(const Point3f &oldDir,
									 const Point3f &normal, Location loc)
{
	Point3f newDir;

	if (loc == Location::Out)
	{
		double cosIN = DotProduct(oldDir, -normal);
		Point3f tmp0 = normal + oldDir/cosIN;
		double cosIN2 = cosIN * cosIN;
		double tmp1 = ri_coef_re + cosIN2 - 1.0;
		tmp1 = (ri_coef_im < FLT_EPSILON) ? tmp1
										  : sqrt(tmp1*tmp1 + ri_coef_im);

		tmp1 = (ri_coef_re + 1.0 - cosIN2 + tmp1)/2.0;
		tmp1 = (tmp1/cosIN2) - Norm(tmp0);
		tmp1 = sqrt(tmp1);
		newDir = (tmp0/tmp1) - normal;
	}
	else // reflection
	{
		double cosIN = DotProduct(oldDir, normal);
		Point3f tmp = oldDir/cosIN - normal;
		newDir = tmp - normal;
	}

	Normalize(newDir);
	return newDir;
}

/** NOTE: Result beams are ordered in inverse direction */
void Tracing::SetPolygonByFacet(int facetId, Polygon &polygon) const
{
	const Polygon &facet = m_facets[facetId];
	int size = facet.size;
	polygon.size = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.arr[i] = facet.arr[size-i];
	}
}

double Tracing::GetIncomingEnergy() const
{
	return m_incommingEnergy;
}

double Tracing::ComputeInternalOpticalPath(const Beam &beam, const vector<int> &tr)
{
#ifdef _DEBUG // DEB
	if (beam.id == /*4194*/11529)
		int ff = 0;
#endif
	Point3f endPoint = /*beam.arr[0]*/beam.Center();
	double path = 0;
	Point3f p1 = endPoint;
	Point3f p2;
	Point3f dir = -beam.light.direction;
	int last = tr.size()-1;

	// first iteration
	dir = ChangeBeamDirection(dir, m_facets[tr[last]].ex_normal, Location::Out);
	ProjectPointToFacet(p1, dir, m_facets[tr[last-1]].in_normal, p2);
	double len = Length(p2 - p1);
#ifdef _DEBUG // DEB
	len *= 1.3116;
#endif

	path += len;
	p1 = p2;

	for (int i = last-2; i >= 0; --i)
	{
		Location loc = beam.GetLocationByLevel(i);
		dir = ChangeBeamDirection(dir, m_facets[tr[i+1]].ex_normal, loc);
		ProjectPointToFacet(p1, dir, m_facets[tr[i]].in_normal, p2);

		if (loc == Location::In)
		{	// sum internal path only
			len = Length(p2 - p1);
#ifdef _DEBUG // DEB
			len *= 1.3116;
#endif
			path += len;
		}

		p1 = p2;
	}
#ifdef _DEBUG // DEB
//	path *= real(m_refrIndex);

	// accounting of far zones
	Point3f pFar1;
	Point3f nFar1 = m_incidentDir;
	nFar1.d_param = FAR_ZONE_DISTANCE;
	ProjectPointToFacet(p1, -nFar1, nFar1, pFar1);

	Point3f pFar2;
	Point3f nFar2 = -beam.light.direction;
	nFar2.d_param = FAR_ZONE_DISTANCE;
	ProjectPointToFacet(endPoint, -nFar2, nFar2, pFar2);

//	double dd2 = Length(p2 - pFar1) /*+ Length(endPoint - pFar2)*/;
	double dd = FAR_ZONE_DISTANCE + DotProductD(p2, nFar1) +
			fabs(DotProductD(endPoint, nFar2) + FAR_ZONE_DISTANCE);
	path += dd;
//cout << beam.opticalPath/*-dd*/ << " " << path;

#endif
#ifdef _DEBUG // DEB
	if (fabs(path - beam.opticalPath) > 7)
		int ff = 0;
#endif
	return path;
}

void Tracing::ProjectPointToFacet(const Point3f &point, const Point3f &direction,
								  const Point3f &facetNormal, Point3f &projection)
{
	double t = DotProduct(point, facetNormal);
	t = t + facetNormal.d_param;
	double dp = DotProduct(direction, facetNormal);
	t = t/dp;
	projection = point - (direction * t);
}
