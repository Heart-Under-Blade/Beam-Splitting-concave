#include "Tracing.h"

#include <float.h>
#include <assert.h>

#include "macro.h"
#include "geometry_lib.h"

#ifdef _DEBUG // DEB
#include <iostream>
using namespace std;
#endif

#define NORM_CEIL	FLT_EPSILON + 1

Tracing::Tracing(Particle *particle, const Point3f &incidentDir, bool isOpticalPath,
				 const Point3f &polarizationBasis, int interReflectionNumber)
{
//	LOG_ASSERT(incidentDir.cx <= NORM_CEIL
//		   && incidentDir.cy <= NORM_CEIL
//		   && incidentDir.cz <= NORM_CEIL
//		   && "Direction of the start beam is not normalized.");

	m_particle = particle;
	m_facets = m_particle->facets;

	m_isOpticalPath = isOpticalPath;
	m_polarizationBasis = polarizationBasis;
	m_interReflectionNumber = interReflectionNumber;
	m_incidentDir = incidentDir;

	m_refrIndex = m_particle->GetRefractionIndex();
	double re = real(m_refrIndex);
	double im = imag(m_refrIndex);
	ri_coef_re = re*re - im*im;
	ri_coef_im = 4*re*re*im;
}

double Tracing::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f normal = m_facets[beam.lastFacetID].ex_normal; // normal of last facet of beam
	double cosFB = DotProduct(normal, beam.direction);
	double e = fabs(cosFB);

	if (e < eps)
	{
		return 0;
	}

	double area = beam.Area();
	double n = Length(normal);
	return (e*area)/n;
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
	SetBeam(outBeam, inBeam, refrDir, inBeam.e, Tv, Th);

	double cosInc2 = (2.0*cosIN);
	SetBeam(inBeam, inBeam, reflDir, inBeam.e,
			cosInc2/Tv0, cosInc2/Th0);
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
	const Point3f &normal = m_facets[facetId].in_normal;
	double cosIN = DotProduct(m_incidentDir, normal);

	if (cosIN < EPS_COS_00) // slopping incidence
	{
		SetSloppingBeamParams_initial(m_incidentDir, cosIN, facetId, inBeam, outBeam);
	}
	else // normal incidence
	{
		inBeam.J.m11 = 2.0/(m_refrIndex + 1.0);
		inBeam.J.m22 = inBeam.J.m11;
		inBeam.e = m_polarizationBasis;
		inBeam.direction = m_incidentDir;

		outBeam.J.m11 = (m_refrIndex - 1.0)/(m_refrIndex + 1.0);
		outBeam.J.m22 = -outBeam.J.m11;
		outBeam.e = m_polarizationBasis;
		outBeam.direction = -m_incidentDir;
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
	beam.e = m_polarizationBasis;
	beam.direction = dir;
	beam.RotatePlane(newBasis);
}

void Tracing::CalcOpticalPath_initial(Beam &inBeam, Beam &outBeam)
{
	Point3f center = inBeam.Center();

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = FAR_ZONE_DISTANCE + DotProduct(m_incidentDir, center);

	outBeam.D = DotProduct(-outBeam.direction, center);
	outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
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

		/// first incident beam
		{
			Beam outBeam;
			TraceFirstBeam(facetId, incidentBeam, outBeam);
			outBuff.push_back(outBeam);
		}

		unsigned int size = tracks.at(i).size();

		try /// internal beams
		{
			for (unsigned int j = 1; j < size; ++j)
			{
				facetId = tracks.at(i).at(j);

				Beam inBeam;
				TraceSecondaryBeams(incidentBeam, facetId, inBeam, outBuff);

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

void Tracing::CalcOpticalPathInternal(double cosIN, const Beam &incidentBeam,
									  Beam &outBeam, Beam &inBeam) const
{
	double Nr = CalcNr(cosIN);
	double coef = (incidentBeam.location == Location::Out) ? 1 : sqrt(Nr);
	Point3f center = outBeam.Center();

	outBeam.D = DotProduct(-outBeam.direction, center);

	double temp = DotProduct(incidentBeam.direction, center);
	outBeam.opticalPath = incidentBeam.opticalPath
			+ coef*fabs(temp + incidentBeam.D);

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = outBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + inBeam.D);
}

bool Tracing::IsTerminalBeam(const Beam &beam)
{
	double j_norm = beam.J.Norm(); // OPT: move to last element of comparison
	return (j_norm < LOW_ENERGY_LEVEL) || (beam.level >= m_interReflectionNumber);
}

double Tracing::CalcNr(const double &cosIN) const
{
	double cosIN_sqr = cosIN*cosIN;
	const double &re = ri_coef_re;
	return (re + sqrt(re*re + ri_coef_im/cosIN_sqr))/2.0;
}

void Tracing::TraceSecondaryBeams(Beam &incidentBeam, int facetID,
								  Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.direction;

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

	if (cosIN >= EPS_COS_00) /// normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);

		outBeam.id = incidentBeam.id;
		outBeam.lastFacetID = facetID;
		outBeam.level = incidentBeam.level + 1;
		SetBeamID(outBeam);
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
			outBeams.push_back(outBeam);
		}
	}
}

void Tracing::SetSloppingIncidenceBeamParams(double cosIN, const Point3f &normal,
											 Beam &incidentBeam,
											 Beam &inBeam, Beam &outBeam,
											 bool &isTrivialIncidence)
{
	const Point3f &incidentDir = incidentBeam.direction;
	double cosIN_sqr = cosIN*cosIN;

	Point3f scatteringNormal;
	CrossProduct(normal, incidentDir, scatteringNormal);
	Normalize(scatteringNormal);
	incidentBeam.RotatePlane(scatteringNormal);

	Point3f r0 = incidentDir/cosIN - normal;

	Point3f reflDir = r0 - normal;
	Normalize(reflDir);

	inBeam.direction = reflDir;
	inBeam.e = scatteringNormal;

	double Nr = CalcNr(cosIN);
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
	const Point3f &dir = incidentBeam.direction;
	complex temp;

	temp = (2.0*m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(outBeam, incidentBeam, dir, incidentBeam.e, temp, temp);

	temp = (1.0 - m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(inBeam, incidentBeam, -dir, incidentBeam.e, temp, -temp);

	if (m_isOpticalPath)
	{
		CalcOpticalPathInternal(cosIN, incidentBeam, inBeam, outBeam);
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
	outBeam.direction = refrDir;
	outBeam.e = inBeam.e;

	complex Tv = (cosIN - tmp1)/Tv0;
	complex Th = (tmp0 - cosRefr)/Th0;
	inBeam.SetJonesMatrix(incidentBeam, Tv, Th);

	if (m_isOpticalPath)
	{
		CalcOpticalPathInternal(Nr, incidentBeam, inBeam, outBeam);
	}
}

void Tracing::SetCompleteReflectionBeamParams(double cosIN, double Nr,
											  const Beam &incidentBeam,
											  Beam &inBeam)
{
	const Point3f &incidentDir = incidentBeam.direction;

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
		inBeam.D = DotProduct(-center, inBeam.direction);

		double temp = DotProduct(incidentDir, center);

		inBeam.opticalPath = incidentBeam.opticalPath
				+ sqrt(Nr)*fabs(temp + incidentBeam.D);
	}
}

void Tracing::SetBeam(Beam &beam, const Beam &other,
					  const Point3f &dir, const Point3f &e,
					  const complex &coef1, const complex &coef2) const
{
	beam.SetJonesMatrix(other, coef1, coef2);
	beam.direction = dir;
	beam.e = e;
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

	bool isProjected = ProjectToFacetPlane(beam, beam.direction, normal,
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
