#include "Tracing.h"

#include <float.h>
#include <assert.h>

#include "macro.h"
#include "geometry_lib.h"

#define NORM_CEIL	FLT_EPSILON + 1

Tracing::Tracing(Particle *particle, const Point3f &startBeamDir, bool isOpticalPath,
				 const Point3f &polarizationBasis, int interReflectionNumber)
{
	LOG_ASSERT(startBeamDir.cx <= NORM_CEIL
		   && startBeamDir.cy <= NORM_CEIL
		   && startBeamDir.cz <= NORM_CEIL
		   && "Direction of the start beam is not normalized.");

	m_particle = particle;
	m_facets = m_particle->facets;

	m_isOpticalPath = isOpticalPath;
	m_polarizationBasis = polarizationBasis;
	m_interReflectionNumber = interReflectionNumber;
	m_initialBeam.direction = startBeamDir;

	m_refrIndex = m_particle->GetRefractionIndex();
	double re = real(m_refrIndex);
	double im = imag(m_refrIndex);
	ri_coef_re = re*re - im*im;
	ri_coef_im = 4*re*re*im;
}

void Tracing::RotateParticle(double beta, double gamma)
{
	m_particle->Rotate(beta, gamma, 0);
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

void Tracing::SetFirstBeamOpticalParams(int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f &startDir = m_initialBeam.direction;
	const Point3f &normal = m_particle->facets[facetId].in_normal;
	double cosIN = DotProduct(startDir, normal);

	if (cosIN < EPS_COS_00) // slopping incidence
	{
		SetSloppingBeamParams_initial(startDir, cosIN, facetId, inBeam, outBeam);
	}
	else // normal incidence
	{
		inBeam.JMatrix.m11 = 2.0/(m_refrIndex + 1.0);
		inBeam.JMatrix.m22 = inBeam.JMatrix.m11;
		inBeam.e = m_polarizationBasis;
		inBeam.direction = startDir;

		outBeam.JMatrix.m11 = (m_refrIndex - 1.0)/(m_refrIndex + 1.0);
		outBeam.JMatrix.m22 = -outBeam.JMatrix.m11;
		outBeam.e = m_polarizationBasis;
		outBeam.direction = -startDir;
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
	Point3f center = CenterOfPolygon(inBeam.polygon);

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = FAR_ZONE_DISTANCE - DotProduct(m_initialBeam.direction, center);

	outBeam.D = DotProduct(-outBeam.direction, center);
	outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.D);
}

void Tracing::SplitExternalBeamByFacet(int facetId, Beam &inBeam, Beam &outBeam)
{
	SetFirstBeamOpticalParams(facetId, inBeam, outBeam);
	SetPolygonByFacet(facetId, inBeam.polygon);
	SetPolygonByFacet(facetId, outBeam.polygon);
}

void Tracing::CalcLigthSurfaceArea(int facetId, const Beam &beam)
{
	const Point3f &startDir = m_initialBeam.direction;
	const Point3f &normal = m_particle->facets[facetId].in_normal;
	double cosIN = DotProduct(startDir, normal);
	m_lightSurfaceArea += AreaOfBeam(beam) * cosIN;
}

// TODO: пофиксить
void Tracing::SplitBeamByParticle(const std::vector<std::vector<int>> &tracks,
								  std::vector<Beam> &outBeams)
{
	for (unsigned int i = 0; i < tracks.size(); ++i)
	{
		int facetId = tracks.at(i).at(0);
		const Point3f &extNormal = m_particle->facets[facetId].ex_normal;

		double cosIN = DotProduct(m_initialBeam.direction, extNormal);

		if (cosIN < EPS_COS_90) /// beam is not incident to this facet
		{
			continue;
		}

		std::vector<Beam> outBuff;
		Beam incidentBeam;

		/// first incident beam
		{
			Beam outBeam;

			SplitExternalBeamByFacet(facetId, incidentBeam, outBeam);

			outBuff.push_back(outBeam);
		}

		unsigned int size = tracks.at(i).size();

		try /// internal beams
		{
			for (unsigned int j = 1; j < size; ++j)
			{
				facetId = tracks.at(i).at(j);

				Beam inBeam;
				SplitInternalBeamByFacet(incidentBeam, facetId, inBeam, outBuff);

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
	double coef = (incidentBeam.location == Location::Outside) ? 1 : sqrt(Nr);
	Point3f center = CenterOfPolygon(outBeam.polygon);

	outBeam.D = DotProduct(-outBeam.direction, center);

	double temp = DotProduct(incidentBeam.direction, center);
	outBeam.opticalPath = incidentBeam.opticalPath
			+ coef*fabs(temp + incidentBeam.D);

	inBeam.D = DotProduct(-inBeam.direction, center);
	inBeam.opticalPath = outBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + inBeam.D);
}

bool Tracing::isTerminalBeam(const Beam &beam)
{
	double j_norm = beam.JMatrix.Norm();
	return (j_norm < LOW_ENERGY_LEVEL) || (beam.level >= m_interReflectionNumber);
}

double Tracing::CalcNr(const double &cosIN) const
{
	double cosIN_sqr = cosIN*cosIN;
	const double &re = ri_coef_re;
	return (re + sqrt(re*re + ri_coef_im/cosIN_sqr))/2.0;
}

void Tracing::SplitInternalBeamByFacet(Beam &incidentBeam, int facetIndex,
									   Beam &inBeam, std::vector<Beam> &outBeams)
{
	Beam outBeam;
	const Point3f &incidentDir = incidentBeam.direction;

	// ext. normal uses in this calculating
	const Point3f &normal = m_particle->facets[facetIndex].ex_normal;
	double cosIN = DotProduct(normal, incidentDir);

	if (cosIN < EPS_COS_90) /// beam is not incident to this facet
	{
		throw std::exception();
	}

	bool isOk = Intersect(facetIndex, incidentBeam, outBeam.polygon);

	if (!isOk)
	{
		throw std::exception();
	}

	inBeam = outBeam;

	if (cosIN >= EPS_COS_00) /// normal incidence
	{
		SetNormalIncidenceBeamParams(cosIN, incidentBeam, inBeam, outBeam);
		outBeams.push_back(outBeam);
	}
	else /// slopping incidence
	{
		bool isTrivialIncidence;
		SetSloppingIncidenceBeamParams(cosIN, normal, incidentBeam,
									   inBeam, outBeam, isTrivialIncidence);
		if (isTrivialIncidence)
		{
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
	const Point3f &incidentDir = incidentBeam.direction;
	complex temp;

	temp = (2.0*m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(outBeam, incidentBeam, incidentDir, incidentBeam.e,
			temp, temp);

	temp = (1.0 - m_refrIndex)/(1.0 + m_refrIndex); // OPT: вынести целиком
	SetBeam(inBeam, incidentBeam, -incidentDir, incidentBeam.e,
			temp, -temp);

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

	MulJMatrix(outBeam, incidentBeam, tmp/Tv0, tmp/Th0);
	outBeam.direction = refrDir;
	outBeam.e = inBeam.e;

	complex Tv = (cosIN - tmp1)/Tv0;
	complex Th = (tmp0 - cosRefr)/Th0;
	MulJMatrix(inBeam, incidentBeam, Tv, Th);

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

	MulJMatrix(inBeam, incidentBeam, Rv, Rh);

	if (m_isOpticalPath)
	{
		Point3f center = CenterOfPolygon(inBeam.polygon);
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
	MulJMatrix(beam, other, coef1, coef2);
	beam.direction = dir;
	beam.e = e;
}

void Tracing::Difference(const Polygon &clip, const Point3f &normal,
						 const Polygon &subject, const Point3f &subjectDir,
						 Polygon *difference, int &resultSize) const
{
	__m128 _clip[MAX_VERTEX_NUM];
	bool isProjected = ProjectToFacetPlane(clip, subjectDir, normal, _clip);

	if (!isProjected)
	{
		difference[resultSize++] = subject;
		return;
	}

	__m128 _clip_normal = _mm_setr_ps(normal.cx, normal.cy, normal.cz, 0.0);
//	__m128 _clip_normal = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
	int clipSize = clip.size;
	__m128 _res_pol[MAX_VERTEX_NUM]; // OPT: заменить на Polygon

	__m128 _subject[MAX_VERTEX_NUM];
	__m128 _buffer[MAX_VERTEX_NUM];

	for (int i = 0; i < subject.size; ++i)
	{
		_subject[i] = _mm_load_ps(subject.arr[i].point);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.size;

	__m128 _first_p, _second_p; // points of projection
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (int i = 0; i < clip.size; ++i)
	{
		int resSize = 0;

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
					__m128 x = computeIntersection_i(_first_p, _second_p, _p1, _p2,
													 _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_res_pol[resSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_buff[bufSize++] = _second_p;
			}
			else
			{
				if (isInFirst)
				{
					__m128 x = computeIntersection_i(_first_p, _second_p, _p1, _p2,
													 _clip_normal, isIntersected);

					if (isIntersected && is_layOnLine_i(x, _first_p, _second_p))
					{
						_res_pol[resSize++] = x;
						_buff[bufSize++] = x;
					}
				}

				_res_pol[resSize++] = _second_p;
			}

			_first_p = _second_p;
			isInFirst = isInSecond;
		}

		if (resSize >= MIN_VERTEX_NUM)
		{
			Polygon resPolygon;
			SetOutputPolygon(_res_pol, resSize, resPolygon);

			if (resPolygon.size >= MIN_VERTEX_NUM)
			{
				difference[resultSize++] = resPolygon;
			}
		}
	}

	if (resultSize == 1)//DEB
		int fff= 0;
}

/// Projection of beam to facet
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

void Tracing::MulJMatrix(Beam &beam1, const Beam &beam2,
						 const complex &coef1, const complex &coef2) const
{
	beam1.JMatrix.m11 = coef1 * beam2.JMatrix.m11;
	beam1.JMatrix.m12 = coef1 * beam2.JMatrix.m12;
	beam1.JMatrix.m21 = coef2 * beam2.JMatrix.m21;
	beam1.JMatrix.m22 = coef2 * beam2.JMatrix.m22;
}

/// NOTE: вершины пучка и грани должны быть ориентированы в одном направлении
bool Tracing::Intersect(int facetId, const Beam &beam, Polygon &intersection) const
{
	__m128 _output_points[MAX_VERTEX_NUM];
	const Point3f &normal = m_particle->facets[facetId].in_normal;

	bool isProjected = ProjectToFacetPlane(beam.polygon, beam.direction,
										   normal, _output_points);
	if (!isProjected)
	{
		return false;
	}

	__m128 _normal_to_facet = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
	__m128 *_output_ptr = _output_points;
	int outputSize = beam.polygon.size;

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	int facetSize = m_particle->facets[facetId].polygon.size;

	__m128 _p1, _p2; // vertices of facet
	__m128 _s_point, _e_point;	// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = m_particle->facets[facetId].polygon.arr[facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = m_particle->facets[facetId].polygon.arr[i];
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
					__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2,
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
				__m128 x = computeIntersection_i(_s_point, _e_point, _p1, _p2,
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

	LOG_ASSERT(tmp1 > 0);
	tmp1 = sqrt(tmp1);

	reflDir = (tmp0/tmp1) - normal;
	Normalize(reflDir);
}

/** NOTE: Result beams are ordered in inverse direction */
void Tracing::SetPolygonByFacet(int facetId, Polygon &polygon) const
{
	const Polygon &facet = m_facets[facetId].polygon;
	int size = facet.size;
	polygon.size = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.arr[i] = facet.arr[size-i];
	}
}

double Tracing::AreaOfBeam(const Beam &beam) const
{
	double square = 0;
	const Point3f &basePoint = beam.polygon.arr[0];
	Point3f p1 = beam.polygon.arr[1] - basePoint;

	for (int i = 2; i < beam.polygon.size; ++i)
	{
		Point3f p2 = beam.polygon.arr[i] - basePoint;
		Point3f res;
		CrossProduct(p1, p2, res);
		square += sqrt(Norm(res));
		p1 = p2;
	}

	return square / 2.0;
}

double Tracing::GetLightSurfaceArea() const
{
	return m_lightSurfaceArea;
}
