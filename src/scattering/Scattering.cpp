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
	: m_particle(particle),
	  m_splitting(isOpticalPath),
	  m_incidentLight(incidentLight),
	  m_nActs(nActs)
{
	m_incidentDir = m_incidentLight->direction;
	m_incidentDir.d_param = m_incidentLight->direction.d_param;

	m_polarBasis = m_incidentLight->polarizationBasis;

	m_splitting.ComputeRiParams(m_particle->GetRefractiveIndex());
}

// REF: перенести в Tracks
IdType Scattering::Scattering::RecomputeTrackId(const IdType &oldId, int facetId)
{
	return (oldId + (facetId + 1)) * (m_particle->nElems + 1);
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
	beam.SetTracingParams(facetId, level, location);
#ifdef _DEBUG // DEB
	beam.dirs.push_back(beam.direction);
#endif
	m_propagatingBeams[m_treeSize++] = beam;
}

void Scattering::SetIncidentBeamOpticalParams(Facet *facet,
											  Beam &inBeam, Beam &outBeam)
{
	const Point3f &normal = facet->in_normal;
	m_splitting.ComputeCosA(m_incidentDir, normal);

	if (!m_splitting.IsNormalIncidence()) // regular incidence
	{
		Beam fakeIncidentBeam;
		fakeIncidentBeam.SetLight(*m_incidentLight);
		const Point3f &facetNormal = facet->in_normal;
		ComputePolarisationParams(-fakeIncidentBeam.direction, facetNormal,
								  fakeIncidentBeam);
		m_splitting.ComputeRegularBeamParamsExternal(facetNormal,
													 fakeIncidentBeam,
													 inBeam, outBeam);
	}
	else // normal incidence
	{
		m_splitting.ComputeNormalBeamParamsExternal(*m_incidentLight,
													inBeam, outBeam);
	}
}

void Scattering::ComputePolarisationParams(const Vector3f &dir,
										   const Vector3f &facetNormal, Beam &beam)
{
	Point3f newBasis = CrossProduct(facetNormal, dir);
	Normalize(newBasis);
	beam.RotateJMatrix(newBasis);
	beam.polarizationBasis = newBasis;
}

void Scattering::SplitLightToBeams(Facet *facet, Beam &inBeam, Beam &outBeam)
{
	SetPolygonByFacet(facet, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
	SetPolygonByFacet(facet, outBeam);
	SetIncidentBeamOpticalParams(facet, inBeam, outBeam);

//	if (m_isOpticalPath)
	{
		Point3f p = inBeam.Center();
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
		inBeam.opticalPath = 0;
		outBeam.opticalPath = 0;
		inBeam.AddOpticalPath(path);
		outBeam.AddOpticalPath(path);
	}
}

Particle *Scattering::GetParticle() const
{
	return m_particle;
}

void Scattering::ComputeFacetEnergy(const Vector3f &facetNormal,
									const Polygon &lightedPolygon)
{
	double cosA = DotProduct(m_incidentDir, facetNormal);
	m_incidentEnergy += lightedPolygon.Area() * cosA;
}

// TODO: пофиксить
void Scattering::ScatterLight(const std::vector<std::vector<int>> &/*tracks*/,
								  std::vector<Beam> &/*outBeams*/)
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

void Scattering::RotateParticle(const Angle &angle)
{
	m_particle->Rotate(angle);
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
	return (beam.nActs >= m_nActs) || (beam.J.Norm() < EPS_BEAM_ENERGY);
}

void Scattering::Difference(const Polygon &subject, const Vector3f &subjNormal,
						 const Polygon &clip, const Vector3f &clipNormal,
						 const Vector3f &clipDir, PolygonArray &difference) const
{
	__m128 _clip[MAX_VERTEX_NUM];
	bool isProjected = ProjectToFacetPlane(clip, clipDir, subjNormal, _clip);

	if (!isProjected)
	{
		difference.Push(subject);
		return;
	}

	__m128 _clip_normal = _mm_setr_ps(clipNormal.cx, clipNormal.cy, clipNormal.cz, 0.0);

	int clipSize = clip.nVertices;
	__m128 _diff_pol[MAX_VERTEX_NUM];

	__m128 _subject[MAX_VERTEX_NUM];
	__m128 _buffer[MAX_VERTEX_NUM];

	for (size_t i = 0; i < subject.nVertices; ++i)
	{
		_subject[i] = _mm_load_ps(subject.arr[i].point);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.nVertices;

	__m128 _first_p, _second_p;
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (size_t i = 0; i < clip.nVertices; ++i)
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

			if (resPolygon.nVertices >= MIN_VERTEX_NUM)
			{
				difference.Push(resPolygon);
			}
		}
	}
}

bool Scattering::ProjectToFacetPlane(const Polygon &polygon, const Vector3f &dir,
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

	for (int i = 0; i < polygon.nVertices; ++i)
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
void Scattering::Intersect(Facet *facet, const Beam &beam, Polygon &intersection) const
{
	__m128 _output_points[MAX_VERTEX_NUM];
	// REF: перенести в случай невыпуклых частиц
	const Point3f &normal = facet->in_normal;
	const Point3f &normal1 = (beam.location == Location::In) ? facet->in_normal
															 : facet->ex_normal;

	bool isProjected = ProjectToFacetPlane(beam, beam.direction, normal1,
										   _output_points);
	if (!isProjected)
	{
		return;
	}

	__m128 _normal_to_facet = _mm_setr_ps(-normal.cx, -normal.cy, -normal.cz, 0.0);
	__m128 *_output_ptr = _output_points;
	int outputSize = beam.nVertices;

	__m128 _buffer[MAX_VERTEX_NUM];
	__m128 *_buffer_ptr = _buffer;
	int bufferSize;

	int facetSize = facet->nVertices;

	__m128 _p1, _p2; // vertices of facet
	__m128 _s_point, _e_point;	// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = facet->arr[facetSize-1];
	_p2 = _mm_load_ps(p2.point);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = facet->arr[i];
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
			polygon.arr[polygon.nVertices++] = p;
		}

		p0 = _output_points[i];
	}
}

/** NOTE: Result beams are ordered in inverse direction */
void Scattering::SetPolygonByFacet(Facet *facet, Polygon &polygon) const
{
	int size = facet->nVertices;
	polygon.nVertices = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.arr[i] = facet->arr[size-i];
	}
}

double Scattering::GetIncedentEnergy() const
{
	return m_incidentEnergy;
}

OpticalPath Scattering::ComputeOpticalPath(const Beam &beam,
										   const Point3f &startPoint,
										   std::vector<int> track)
{
	OpticalPath path;

	Point3f dir = -beam.direction; // back direction
	Location loc = Location::Out;
	Location nextLoc;

	Point3f p1 = startPoint;
	Point3f p2;

	// back tracing
	for (int i = track.size()-1; i > 0; --i)
	{
		nextLoc = beam.GetLocationByActNumber(i-1);

		const Point3f &exNormal = m_particle->GetActualFacet(track[i])->ex_normal;
		dir = m_splitting.ChangeBeamDirection(dir, exNormal, loc, nextLoc);

		const Point3f &inNormal = m_particle->GetActualFacet(track[i-1])->in_normal;
		p2 = ProjectPointToPlane(p1, dir, inNormal);
		double len = Length(p2 - p1);

		if (nextLoc == Location::In)
		{	// add internal path only
#ifdef _DEBUG // DEB
			m_splitting.ComputeCosA(dir, exNormal);
			double reRi = m_splitting.ComputeEffectiveReRi();
			len *= sqrt(reRi);
#endif
			path.internal += len;
		}
		else
		{
			path.external += len;
		}

		p1 = p2;
		loc = nextLoc;
	}

#ifdef _DEBUG // DEB
//	path *= real(m_splitting.GetRi());
	Point3f nFar1 = m_incidentDir;
	Point3f nFar2 = -beam.direction;
	double dd1 = m_splitting.FAR_ZONE_DISTANCE + DotProductD(p2, nFar1);
	double dd2 = fabs(DotProductD(beam.Center(), nFar2) + m_splitting.FAR_ZONE_DISTANCE);

	path.external += dd1;
	path.external += dd2;

//	if (fabs(path.GetTotal() - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}
