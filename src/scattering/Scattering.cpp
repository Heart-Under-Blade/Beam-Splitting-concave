#include "Scattering.h"

#include <float.h>
#include <assert.h>

#include "macro.h"
#include "geometry_lib.h"

#ifdef _DEBUG // DEB
#include <iostream>
#endif

#define NORM_CEIL	FLT_EPSILON + 1

using namespace std;

Scattering::Scattering(Particle *particle, Light *incidentLight, bool isOpticalPath,
					   int nActs)
	: m_incidentLight(incidentLight),
	  m_particle(particle)
{
	m_facets = m_particle->facets;

	m_isOpticalPath = isOpticalPath;
	m_nActs = nActs;

	m_incidentDir = m_incidentLight->direction;
	m_incidentDir.d_param = m_incidentLight->direction.d_param;

	m_polarBasis = m_incidentLight->polarizationBasis;

	m_ri = m_particle->GetRefractiveIndex();
	double re = real(m_ri);
	double im = imag(m_ri);
	m_cRiRe = re*re - im*im;
	m_cRiRe2 = m_cRiRe * m_cRiRe;
	m_cRiIm = 4*re*re*im;
}

void Scattering::SetRegularBeamParamsExternal(const Point3f &beamDir, double cosA,
											  int facetId, Beam &inBeam, Beam &outBeam)
{
	const Point3f &facetNormal = m_facets[facetId].in_normal;

	RotatePolarisationPlane(beamDir, facetNormal, inBeam);

	Point3f refrDir, reflDir;
	SplitDirection(beamDir, cosA, -facetNormal, reflDir, refrDir);

	double cosB = DotProduct(facetNormal, reflDir);

	complex Tv00 = m_ri * cosA;
	complex Th00 = m_ri * cosB;

	complex Tv0 = Tv00 + cosB;
	complex Th0 = Th00 + cosA;

	complex Tv = (Tv00 - cosB)/Tv0;
	complex Th = (cosA - Th00)/Th0;

	outBeam.SetJonesMatrix(inBeam, Tv, Th);
	outBeam.SetLight(refrDir, inBeam.polarizationBasis);

	double cosA2 = 2.0*cosA;

	inBeam.SetJonesMatrix(inBeam, cosA2/Tv0, cosA2/Th0);
	inBeam.direction = reflDir;
}

void Scattering::ComputeBeamId(Beam &beam)
{
	beam.trackId += (beam.lastFacetId + 1);
	beam.trackId *= (m_particle->nFacets + 1);
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
	beam.SetTracingParams(facetId, level, location);
	PushBeamToTree(beam);
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level)
{
	beam.lastFacetId = facetId;
	beam.level = level;
	PushBeamToTree(beam);
}

void Scattering::PushBeamToTree(Beam &beam)
{
	ComputeBeamId(beam);
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
		inBeam = (*m_incidentLight);
		inBeam.J.m11 = 2.0/(m_ri + 1.0);
		inBeam.J.m22 = inBeam.J.m11;

		outBeam.J.m11 = (m_ri - 1.0)/(m_ri + 1.0);
		outBeam.J.m22 = -outBeam.J.m11;
		outBeam.direction = -dir;
		outBeam.polarizationBasis = m_incidentLight->polarizationBasis;
	}

	if (m_isOpticalPath)
	{
		CalcOpticalPathForLight(inBeam, outBeam);
	}
}

void Scattering::RotatePolarisationPlane(const Point3f &dir,
										 const Point3f &facetNormal,
										 Beam &beam)
{
	Point3f newBasis;
	CrossProduct(facetNormal, -dir, newBasis);
	Normalize(newBasis);
	beam.direction = dir;
	beam.polarizationBasis = m_incidentLight->polarizationBasis;
	beam.RotatePlane(newBasis);
}

void Scattering::CalcOpticalPathForLight(Beam &inBeam, Beam &outBeam)
{
	inBeam.frontPosition = DotProduct(-inBeam.direction, inBeam.arr[0]);
	inBeam.opticalPath = FAR_ZONE_DISTANCE + DotProduct(m_incidentDir, inBeam.arr[0]);

	outBeam.frontPosition = DotProduct(-outBeam.direction, inBeam.arr[0]);
	outBeam.opticalPath = inBeam.opticalPath + fabs(FAR_ZONE_DISTANCE + outBeam.frontPosition);
}

void Scattering::SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam)
{
	SetPolygonByFacet(facetId, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
	SetPolygonByFacet(facetId, outBeam);
	SetBeamOpticalParams(facetId, inBeam, outBeam);
}

void Scattering::ComputeFacetEnergy(int facetId, const Polygon &lightedPolygon)
{
	const Point3f &normal = m_facets[facetId].in_normal;
	double cosIN = DotProduct(m_incidentDir, normal);
	m_incommingEnergy += lightedPolygon.Area() * cosIN;
}

// TODO: пофиксить
void Scattering::ScatterLight(double beta, double gamma, const std::vector<std::vector<int>> &tracks,
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

//		// first incident beam
//		{
//			Beam outBeam;
//			SplitLightToBeams(facetId, incidentBeam, outBeam);
//			outBuff.push_back(outBeam);
//		}

//		unsigned int size = tracks.at(i).size();

//		try // internal beams
//		{
//			for (unsigned int j = 1; j < size; ++j)
//			{
//				facetId = tracks.at(i).at(j);

//				Beam inBeam;
//				SplitSecondaryBeams(incidentBeam, facetId, inBeam, outBuff);

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

double Scattering::ComputeOpticalPath(const Beam &beam, double cosA,
									  const Point3f &facetPoint) const
{
	double tmp = DotProduct(beam.direction, facetPoint);
	double path = fabs(tmp + beam.frontPosition); // refractive index of external environment = 1

	if (beam.location == Location::In)
	{
		path *= sqrt(ComputeReRI(cosA));
	}

	path += beam.opticalPath;
}

void Scattering::ComputeOpticalParams(double cosA, const Beam &incidentBeam,
									  Beam &inBeam, Beam &outBeam) const
{
	Point3f &facetPoint = inBeam.arr[0];
	double path = ComputeOpticalPath(incidentBeam, cosA, facetPoint);

	inBeam.frontPosition = DotProduct(-inBeam.direction, facetPoint);
	inBeam.opticalPath = path;

	outBeam.frontPosition = DotProduct(-outBeam.direction, facetPoint);
	outBeam.opticalPath = path;
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
	return (beam.level >= m_nActs) || (beam.J.Norm() < EPS_BEAM_ENERGY);
}

double Scattering::ComputeReRI(const double &cosA) const
{
	return (m_cRiRe + sqrt(m_cRiRe2 + m_cRiIm/(cosA*cosA)))/2.0;
}

void Scattering::SetRegularIncidenceBeamParams(double cosA, const Point3f &normal,
											   Beam &incidentBeam,
											   Beam &inBeam, Beam &outBeam,
											   bool &isRegular)
{
	const Point3f &dir = incidentBeam.direction;

	Point3f scatteringNormal;
	CrossProduct(normal, dir, scatteringNormal);
	Normalize(scatteringNormal);
	incidentBeam.RotatePlane(scatteringNormal);

	Point3f r0 = dir/cosA - normal;

	Point3f reflDir = r0 - normal;
	Normalize(reflDir);

	inBeam.SetLight(reflDir, scatteringNormal);

	double reRI = ComputeReRI(cosA);
	double s = 1.0/(reRI*cosA*cosA) - Norm(r0);

	if (s > DBL_EPSILON) // regular incidence
	{
		SetRegularBeamParams(cosA, normal, r0, s, incidentBeam,
							 inBeam, outBeam);
		isRegular = true;
	}
	else // complete internal reflection incidence
	{
		SetCRBeamParams(cosA, reRI, incidentBeam, inBeam);
		isRegular = false;
	}
}

void Scattering::SetNormalIncidenceBeamParams(double cosA, const Beam &incidentBeam,
											  Beam &inBeam, Beam &outBeam)
{
	const Point3f &dir = incidentBeam.direction;
	complex temp;

	temp = (2.0 * m_ri)/(1.0 + m_ri); // OPT: вынести целиком
	outBeam.SetJonesMatrix(incidentBeam, temp, temp);
	outBeam.SetLight(dir, incidentBeam.polarizationBasis);

	temp = (1.0 - m_ri)/(1.0 + m_ri); // OPT: вынести целиком
	inBeam.SetJonesMatrix(incidentBeam, temp, -temp);
	inBeam.SetLight(-dir, incidentBeam.polarizationBasis);

	if (m_isOpticalPath)
	{
		ComputeOpticalParams(cosA, incidentBeam, inBeam, outBeam);
	}
}

void Scattering::SetRegularBeamParams(double cosA, const Point3f &normal,
									  Point3f r0, double s,
									  const Beam &incidentBeam,
									  Beam &inBeam, Beam &outBeam)
{
	Point3f refrDir = r0/sqrt(s) + normal;
	Normalize(refrDir);

	double cosG = DotProduct(normal, refrDir);

	complex tmp0 = m_ri * cosA;
	complex tmp1 = m_ri * cosG;
	complex tmp = 2.0*tmp0;
	complex Tv0 = tmp1 + cosA;
	complex Th0 = tmp0 + cosG;

	outBeam.SetJonesMatrix(incidentBeam, tmp/Tv0, tmp/Th0);
	outBeam = Light{refrDir, inBeam.polarizationBasis};

	complex Tv = (cosA - tmp1)/Tv0;
	complex Th = (tmp0 - cosG)/Th0;
	inBeam.SetJonesMatrix(incidentBeam, Tv, Th);

	if (m_isOpticalPath)
	{
		ComputeOpticalParams(cosA, incidentBeam, inBeam, outBeam);
	}
}

void Scattering::SetCRBeamParams(double cosA, double reRi,
								 const Beam &incidentBeam, Beam &inBeam)
{
	const Point3f &incidentDir = incidentBeam.direction;

	complex tmp0 = m_ri * cosA;
	const double bf = reRi*(1.0 - cosA*cosA) - 1.0;
	double im = (bf > 0) ? sqrt(bf) : 0;

	const complex sq(0, im);
	complex tmp = m_ri * sq;
	complex Rv = (cosA - tmp)/(tmp + cosA);
	complex Rh = (tmp0 - sq)/(tmp0 + sq);

	inBeam.SetJonesMatrix(incidentBeam, Rv, Rh);

	if (m_isOpticalPath)
	{
		Point3f center = inBeam.Center();
		inBeam.frontPosition = DotProduct(-center, inBeam.direction);

		double temp = DotProduct(incidentDir, center);

		inBeam.opticalPath = incidentBeam.opticalPath
				+ sqrt(reRi)*fabs(temp + incidentBeam.frontPosition);
	}
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

void Scattering::SplitDirection(const Point3f &incidentDir, double cosA,
								const Point3f &normal,
								Point3f &reflDir, Point3f &refrDir) const
{
	Point3f tmp0 = normal + incidentDir/cosA;

	refrDir = normal + tmp0;
	Normalize(refrDir);

	double cosA2 = cosA * cosA;
	double tmp1 = m_cRiRe + cosA2 - 1.0;

	if (m_cRiIm > FLT_EPSILON)
	{
		tmp1 = sqrt(tmp1*tmp1 + m_cRiIm);
	}

	tmp1 = (m_cRiRe + 1.0 - cosA2 + tmp1)/2.0;
	tmp1 = (tmp1/cosA2) - Norm(tmp0);

//	LOG_ASSERT(tmp1 > 0);
	tmp1 = sqrt(tmp1);

	reflDir = (tmp0/tmp1) - normal;
	Normalize(reflDir);
}

Point3f Scattering::ChangeBeamDirection(const Vector3f &oldDir,
										const Vector3f &normal, Location loc)
{
	Point3f newDir;

	if (loc == Location::Out) // refraction
	{
		double cosA = DotProduct(oldDir, -normal);
		double cosA2 = cosA * cosA;

		Vector3f tmp0 = normal + oldDir/cosA;
		double tmp1 = m_cRiRe + cosA2 - 1.0;

		if (m_cRiIm > FLT_EPSILON)
		{
			tmp1 = sqrt(tmp1*tmp1 + m_cRiIm);
		}

		tmp1 = (m_cRiRe + 1.0 - cosA2 + tmp1)/2.0;
		tmp1 = (tmp1/cosA2) - Norm(tmp0);
		tmp1 = sqrt(tmp1);
		newDir = (tmp0/tmp1) - normal;
	}
	else // reflection
	{
		double cosA = DotProduct(oldDir, normal);
		Point3f tmp = oldDir/cosA - normal;
		newDir = tmp - normal;
	}

	Normalize(newDir);
	return newDir;
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

double Scattering::ComputeInternalOpticalPath(const Beam &beam,
											  const vector<int> &track)
{
	double path = 0;
	Point3f dir = -beam.direction; // back direction
	Location loc = Location::Out;

	Point3f p1 = beam.arr[0];
	Point3f p2;

	for (int lvl = track.size()-1; lvl > 0; --lvl)
	{
		Point3f &extexnalNormal = m_facets[track[lvl]].ex_normal;
		dir = ChangeBeamDirection(dir, extexnalNormal, loc);

		Point3f &intexnalNormal = m_facets[track[lvl-1]].in_normal;
		p2 = ProjectPointToPlane(p1, dir, intexnalNormal);

		loc = beam.GetLocationByLevel(lvl-1);

		if (loc == Location::In)
		{	// add internal path only
			path += Length(p1 - p2);
		}

		p1 = p2;
	}

	return path;
}
