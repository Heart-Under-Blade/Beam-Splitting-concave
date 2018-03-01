#include "Scattering.h"

#include <float.h>
#include <assert.h>

#include "macro.h"
#include "geometry_lib.h"

#ifdef _DEBUG // DEB
#include <iostream>
using namespace std;
#endif

#define NORM_CEIL	FLT_EPSILON + 1

Scattering::Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
				 int interReflectionNumber)
	: m_incidentLight(incidentLight),
	  m_particle(particle)
{
	m_facets = m_particle->facets;

	m_isOpticalPath = isOpticalPath;
	m_numberOfActs = interReflectionNumber;

	m_incidentDir = m_incidentLight->direction;
	m_incidentDir.d_param = m_incidentLight->direction.d_param;

	m_polarBasis = m_incidentLight->polarizationBasis;

	m_refrIndex = m_particle->GetRefractiveIndex();
	double re = real(m_refrIndex);
	double im = imag(m_refrIndex);
	ri_coef_re = re*re - im*im;
	ri_coef_im = 4*re*re*im;
}

void Scattering::SetRegularBeamParamsExternal(const Point3f &beamDir, double cosIN,
											int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f &facetNormal = m_particle->facets[facetId].in_normal;

	RotatePolarisationPlane(beamDir, facetNormal, inBeam);

	Point3f refrDir, reflDir;
	SplitDirection(beamDir, cosIN, -facetNormal, reflDir, refrDir);

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
}

void Scattering::SetBeamID(Beam &beam)
{
#ifdef _TRACK_ALLOW
	beam.id += (beam.lastFacetID + 1);
	beam.id *= (m_particle->facetNum + 1);
	//	AddToTrack(beam, facetId);
#endif
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
	beam.SetTracingParams(facetId, level, location);
	PushBeamToTree(beam);
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level)
{
	beam.lastFacetID = facetId;
	beam.level = level;
	PushBeamToTree(beam);
}

void Scattering::PushBeamToTree(Beam &beam)
{
	SetBeamID(beam);
	m_beamTree[m_treeSize++] = beam;
}

void Scattering::SetBeamOpticalParams(int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f dir = m_incidentLight->direction;
	const Point3f &normal = m_facets[facetId].in_normal;
	double cosIN = DotProduct(dir, normal);

	if (cosIN < EPS_COS_00) // regular incidence
	{
		SetRegularBeamParamsExternal(dir, cosIN, facetId, inBeam, outBeam);
	}
	else // normal incidence
	{
		inBeam.J.m11 = 2.0/(m_refrIndex + 1.0);
		inBeam.J.m22 = inBeam.J.m11;
		inBeam.light = (*m_incidentLight);

		outBeam.J.m11 = (m_refrIndex - 1.0)/(m_refrIndex + 1.0);
		outBeam.J.m22 = -outBeam.J.m11;
		outBeam.light = Light{-dir, m_incidentLight->polarizationBasis};
	}

	if (m_isOpticalPath)
	{
		CalcOpticalPathForLight(inBeam, outBeam);
	}
}

void Scattering::RotatePolarisationPlane(const Point3f &dir, const Point3f &facetNormal,
									  Beam &beam)
{
	Point3f newBasis;
	CrossProduct(facetNormal, -dir, newBasis);
	Normalize(newBasis);
	beam.light = Light{dir, m_incidentLight->polarizationBasis};
	beam.RotatePlane(newBasis);
}

void Scattering::CalcOpticalPathForLight(Beam &inBeam, Beam &outBeam)
{
	Point3f center = inBeam.Center();

	inBeam.D = DotProduct(-inBeam.light.direction, center);
	inBeam.opticalPath = FAR_ZONE_DISTANCE + DotProduct(m_incidentDir, center);

	outBeam.D = DotProduct(-outBeam.light.direction, center);
	outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
}

void Scattering::SplitLightIntoBeams(int facetId, Beam &inBeam, Beam &outBeam)
{
	SetPolygonByFacet(facetId, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
	SetPolygonByFacet(facetId, outBeam);
	SetBeamOpticalParams(facetId, inBeam, outBeam);
}

void Scattering::CalcFacetEnergy(int facetID, const Polygon &lightedPolygon)
{
	const Point3f &normal = m_facets[facetID].in_normal;
	double cosIN = DotProduct(m_incidentDir, normal);
	m_incommingEnergy += lightedPolygon.Area() * cosIN;
}

// TODO: пофиксить
void Scattering::ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
								  std::vector<Beam> &outBeams)
{
	m_particle->Rotate(beta, gamma, 0);

	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		int facetId = tracks.at(i).at(0);
		const Point3f &extNormal = m_facets[facetId].ex_normal;

		double cosIN = DotProduct(m_incidentDir, extNormal);

		if (cosIN < EPS_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		std::vector<Beam> outBuff;
		Beam incidentBeam;

		// first incident beam
		{
			Beam outBeam;
			SplitLightIntoBeams(facetId, incidentBeam, outBeam);
			outBuff.push_back(outBeam);
		}

		unsigned int size = tracks.at(i).size();

		try // internal beams
		{
			for (unsigned int j = 1; j < size; ++j)
			{
				facetId = tracks.at(i).at(j);

				Beam inBeam;
				SplitSecondaryBeams(incidentBeam, facetId, inBeam, outBuff);

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

void Scattering::CalcOpticalPath(double cosIN, const Beam &incidentBeam,
								 Beam &inBeam, Beam &outBeam) const
{
	Point3f center = inBeam.Center();

	// refractive index of external environment = 1
	double OP = fabs(DotProduct(incidentBeam.light.direction, center) + incidentBeam.D);

	if (incidentBeam.location == Location::In)
	{
		OP *= sqrt(CalcReRI(cosIN));
		inBeam.internalOpticalPath = outBeam.internalOpticalPath = OP;
	}

	inBeam.D = DotProduct(-inBeam.light.direction, center);
	inBeam.opticalPath = incidentBeam.opticalPath + OP;

	outBeam.D = DotProduct(-outBeam.light.direction, center);
	outBeam.opticalPath = inBeam.opticalPath;
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
	double j_norm = beam.J.Norm(); // OPT: move to last element of comparison
	return (j_norm < LOW_ENERGY_LEVEL) || (beam.level >= m_numberOfActs);
}

double Scattering::CalcReRI(const double &cosIN) const
{
	const double &re = ri_coef_re;
	return (re + sqrt(re*re + ri_coef_im/(cosIN*cosIN)))/2.0;
}

void Scattering::SplitSecondaryBeams(Beam &incidentBeam, int facetID,
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

	if (cosIN >= EPS_COS_00) // normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);

		outBeam.id = incidentBeam.id;
		outBeam.lastFacetID = facetID;
		outBeam.level = incidentBeam.level + 1;
		SetBeamID(outBeam);
		outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.D); // добираем оптический путь
		outBeams.push_back(outBeam);
	}
	else // regular incidence
	{
		bool isTrivialIncidence;
		SetRegularBeamParams(cosIN, normal, incidentBeam,
									   inBeam, outBeam, isTrivialIncidence);
		if (isTrivialIncidence)
		{
			outBeam.id = incidentBeam.id;
			outBeam.lastFacetID = facetID;
			outBeam.level = incidentBeam.level + 1;
			SetBeamID(outBeam);
			outBeam.opticalPath += fabs(FAR_ZONE_DISTANCE + outBeam.D); // добираем оптический путь
			outBeams.push_back(outBeam);
		}
	}
}

void Scattering::SetRegularBeamParams(double cosIN, const Point3f &normal,
											 Beam &incidentBeam,
											 Beam &inBeam, Beam &outBeam,
											 bool &isRegular)
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

	double reRI = CalcReRI(cosIN);
	double s = 1.0/(reRI*cosIN_sqr) - Norm(r0);

	if (s > DBL_EPSILON) // regular incidence
	{
		SetRegularBeamParams(cosIN, reRI, normal, r0, s, incidentBeam,
							 inBeam, outBeam);
		isRegular = true;
	}
	else // complete internal reflection incidence
	{
		SetCRBeamParams(cosIN, reRI, incidentBeam, inBeam);
		isRegular = false;
	}
}

void Scattering::SetNormalIncidenceBeamParams(double cosIN, const Beam &incidentBeam,
										   Beam &inBeam, Beam &outBeam)
{
	const Point3f &dir = incidentBeam.light.direction;
	complex temp;

	temp = (2.0*m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(outBeam, incidentBeam, dir, temp, temp);

	temp = (1.0 - m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(inBeam, incidentBeam, -dir, temp, -temp);

	if (m_isOpticalPath)
	{
		CalcOpticalPath(cosIN, incidentBeam, inBeam, outBeam);
	}
}

void Scattering::SetRegularBeamParams(double cosIN, double reRI,
									  const Point3f &normal, Point3f r0,
									  double s, const Beam &incidentBeam,
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

	complex Tv = (cosIN - tmp1)/Tv0;
	complex Th = (tmp0 - cosRefr)/Th0;
	inBeam.SetJonesMatrix(incidentBeam, Tv, Th);

	if (m_isOpticalPath)
	{
		CalcOpticalPath(cosIN, incidentBeam, inBeam, outBeam);
	}
}

void Scattering::SetCRBeamParams(double cosIN, double Nr,
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
		Point3f center = inBeam.Center();
		inBeam.D = DotProduct(-center, inBeam.light.direction);

		double temp = DotProduct(incidentDir, center);

		inBeam.opticalPath = incidentBeam.opticalPath
				+ sqrt(Nr)*fabs(temp + incidentBeam.D);
	}
}

void Scattering::SetBeam(Beam &beam, const Beam &other, const Point3f &dir,
					  const complex &coef1, const complex &coef2) const
{
	beam.SetJonesMatrix(other, coef1, coef2);
	beam.light = Light{dir, other.light.polarizationBasis};
}

void Scattering::Difference(const Polygon &subject, const Point3f &subjNormal,
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

bool Scattering::ProjectToFacetPlane(const Polygon &polygon, const Point3f &dir,
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
bool Scattering::Intersect(int facetID, const Beam &beam, Polygon &intersection) const
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

void Scattering::SetOutputPolygon(__m128 *_output_points, int outputSize,
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

void Scattering::SplitDirection(const Point3f &incidentDir, double cosIN,
								  const Point3f &normal,
								  Point3f &reflDir, Point3f &refrDir) const
{
	Point3f tmp0 = normal + incidentDir/cosIN;

	refrDir = normal + tmp0;
	Normalize(refrDir);

	double cosI_sqr = cosIN * cosIN;
	double tmp1 = ri_coef_re + cosI_sqr - 1.0;

	tmp1 = (ri_coef_im - FLT_EPSILON < 0) ? tmp1
										  : sqrt(tmp1*tmp1 + ri_coef_im);

	tmp1 = (ri_coef_re + 1.0 - cosI_sqr + tmp1)/2.0;
	tmp1 = (tmp1/cosI_sqr) - Norm(tmp0);

//	LOG_ASSERT(tmp1 > 0);
	tmp1 = sqrt(tmp1);

	reflDir = (tmp0/tmp1) - normal;
	Normalize(reflDir);
}

/** NOTE: Result beams are ordered in inverse direction */
void Scattering::SetPolygonByFacet(int facetId, Polygon &polygon) const
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

double Scattering::GetIncomingEnergy() const
{
	return m_incommingEnergy;
}

//double Tracing::CrossSection(const Point3f &beamDir) const
//{
//	double cs = 0.0;

//	for (int i = 0; i < m_particle->facetNum; ++i)
//	{
//		const Point3f n = m_facets[i].Normal();
//		double csa = DotProduct(beamDir, n);

//		if (csa < EPS_COS_90)
//		{
//			continue;
//		}

//		cs += m_facets[i].Area()*csa;
//	}

//	return cs;
//}
