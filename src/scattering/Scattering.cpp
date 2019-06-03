#include "Scattering.h"

#include <float.h>
#include <assert.h>
#include <algorithm>

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
	m_facets = m_particle->facets;

	m_incidentDir = m_incidentLight->direction;
	m_incidentDir.d_param = m_incidentLight->direction.d_param;

	m_polarBasis = m_incidentLight->polarizationBasis;

	m_splitting.ComputeRiParams(m_particle->GetRefractiveIndex());
}

IdType Scattering::Scattering::RecomputeTrackId(const IdType &oldId, int facetId)
{
	return (oldId + (facetId + 1)) * (m_particle->nFacets + 1);
}

void Scattering::PushBeamToTree(Beam &beam, int facetId, int level, Location location)
{
	beam.SetTracingParams(facetId, level, location);
#ifdef _DEBUG // DEB
	beam.dirs.push_back(beam.direction);
#endif
	m_beamTree[m_treeSize++] = beam;
}

void Scattering::SetIncidentBeamOpticalParams(unsigned facetId,
											  Beam &inBeam, Beam &outBeam)
{
	const Point3f &normal = m_facets[facetId].in_normal;
	m_splitting.ComputeCosA(m_incidentDir, normal);

	if (!m_splitting.IsNormalIncidence()) // regular incidence
	{
		Beam fakeIncidentBeam;
		fakeIncidentBeam.SetLight(*m_incidentLight);
		const Point3f &facetNormal = m_facets[facetId].in_normal;
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

void Scattering::SplitLightToBeams(int facetId, Beam &inBeam, Beam &outBeam)
{
	SetPolygonByFacet(facetId, inBeam); // REF: try to exchange this to inBeam = m_facets[facetId]
	SetPolygonByFacet(facetId, outBeam);
	SetIncidentBeamOpticalParams(facetId, inBeam, outBeam);

//	if (m_isOpticalPath)
	{
		Point3f p = inBeam.Center();
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
		inBeam.opticalPath = 0;
		outBeam.opticalPath = 0;
		inBeam.AddOpticalPath(path);
		outBeam.AddOpticalPath(path);
		outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam);
	}
}

void Scattering::ComputeFacetEnergy(int facetId, const Polygon &lightedPolygon)
{
	const Point3f &normal = m_facets[facetId].in_normal;
	double cosA = DotProduct(m_incidentDir, normal);
	m_incidentEnergy += lightedPolygon.Area() * cosA;
}

// TODO: пофиксить
void Scattering::ScatterLight(double /*beta*/, double /*gamma*/, const std::vector<std::vector<int>> &/*tracks*/,
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

void Scattering::OrderVertices2f(std::vector<Point2f> &vertices,
								 Polygon &orderedPolygon)
{
	// define base point
	int baseIndex = 0;
	double minX = vertices[baseIndex].x;

	for (int i = 1; i < vertices.size(); ++i)
	{
		if (vertices[i].x < minX)
		{
			baseIndex = i;
			minX = vertices[i].x;
		}
	}

	std::swap(vertices[baseIndex], vertices.back());
	Point2f base = vertices.back();
//	orderedPolygon.AddVertex(Point3f(base.x, base.y, 0));

	int iBase = 0;

	for (int i = 0; iBase != vertices.size()-1 && i < vertices.size(); ++i)
	{
		iBase = i;
		Point2f vBase = vertices[iBase] - base;

		for (int j = i + 1; j <= vertices.size(); ++j)
		{
			int iNext = (j == vertices.size()) ? 0 : j;
			Point2f vNext = vertices[iNext] - base;

			if (vBase.CrossProduct(vNext) > 0)
			{
				iBase = iNext;
				vBase = vNext;
			}
		}

		std::swap(vertices[iBase], vertices[i]);
		base = vertices[i];
		orderedPolygon.AddVertex(Point3f(base.x, base.y, 10000));
	}
}

void Scattering::ProjectParticleToXY(std::vector<Point2f> &projected)
{
	Point3f n(0, 0, 1, 10000);
	n.d_param = m_particle->MaximalDimention();

	for (int i = 0; i < m_particle->nFacets; i++)
	{
		auto &f = m_particle->facets[i];

		if (DotProduct(f.in_normal, -n) < EPS_COS_90)
		{
			for (int j = 0; j < f.nVertices; j++)
			{
//				auto p = ProjectPointToPlane(f.arr[j], -m_incidentDir, n);
//				projected.push_back(Point2f(-p.coordinates[0], -p.coordinates[1]));
				double tmp = (n.d_param - DotProduct(n, f.arr[j]));
				auto p = f.arr[j] + n*tmp;
				projected.push_back(Point2f(p.coordinates[0], p.coordinates[1]));
			}
		}
	}
}

void Scattering::RemoveDublicatedVertices2f(const std::vector<Point2f> &projected,
										  std::vector<Point2f> &cleared)
{
	for (int i = 0; i < projected.size(); ++i)
	{
		bool isUnique = true;

		for (int j = i + 1; j < projected.size(); ++j)
		{
			if (projected[i].IsEqualTo(projected[j], 0.0001))
			{
				isUnique = false;
			}
		}

		if (isUnique)
		{
			cleared.push_back(projected[i]);
		}
	}
}

void Scattering::FormShadowBeam(std::vector<Beam> &scaterredBeams)
{
	std::vector<Point2f> projected;
	ProjectParticleToXY(projected);

	std::vector<Point2f> projectedCleared;
	RemoveDublicatedVertices2f(projected, projectedCleared);

	Beam shadowBeam;
	OrderVertices2f(projectedCleared, shadowBeam);

	Matrix2x2c jones;
	jones.m11 = -jones.m11;
	jones.m22 = -jones.m22;
	shadowBeam.SetMatrix(jones);

	shadowBeam.direction = m_incidentLight->direction;
	shadowBeam.polarizationBasis = m_incidentLight->polarizationBasis;
	shadowBeam.opticalPath = 20000;
	shadowBeam.lastFacetId = INT_MAX;
	scaterredBeams.push_back(shadowBeam);
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

	for (int i = 0; i < subject.nVertices; ++i)
	{
		_subject[i] = _mm_load_ps(subject.arr[i].coordinates);
	}

	__m128 *_subj = _buffer;
	__m128 *_buff = _subject;
	int bufSize = subject.nVertices;

	__m128 _first_p, _second_p;
	bool isInFirst, isInSecond;

	__m128 _p2 = _clip[clipSize-1];

	for (int i = 0; i < clip.nVertices; ++i)
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
void Scattering::Intersect(int facetID, const Beam &beam, Polygon &intersection) const
{
	__m128 _output_points[MAX_VERTEX_NUM];
	// REF: перенести в случай невыпуклых частиц
	const Point3f &normal = m_facets[facetID].in_normal;

	const Point3f &normal1 = (beam.location == Location::In) ? m_facets[facetID].in_normal
															: m_facets[facetID].ex_normal;
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

	int facetSize = m_particle->facets[facetID].nVertices;

	__m128 _p1, _p2; // vertices of facet
	__m128 _s_point, _e_point;	// points of projection
	bool isInsideE, isInsideS;

	Point3f p2 = m_particle->facets[facetID].arr[facetSize-1];
	_p2 = _mm_load_ps(p2.coordinates);

	for (int i = 0; i < facetSize; ++i)
	{
		_p1 = _p2;
		p2 = m_particle->facets[facetID].arr[i];
		_p2 = _mm_load_ps(p2.coordinates);

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
void Scattering::SetPolygonByFacet(int facetId, Polygon &polygon) const
{
	const Polygon &facet = m_facets[facetId];
	int size = facet.nVertices;
	polygon.nVertices = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.arr[i] = facet.arr[size-i];
	}
}

double Scattering::GetIncedentEnergy() const
{
	return m_incidentEnergy;
}

double Scattering::ComputeInternalOpticalPath(const Beam &beam,
											  const Point3f sourcePoint,
											  const vector<int> &track)
{
	double path1 = 0;
	double path = 0;
	Point3f dir = -beam.direction; // back direction
	Location loc = Location::Out;
	Location nextLoc;

	Point3f p1 = sourcePoint;
	Point3f p2;

	for (int i = track.size()-1; i > 0; --i)
	{
		nextLoc = beam.GetLocationByActNumber(i-1);

		Point3f &exNormal = m_facets[track[i]].ex_normal;
		dir = m_splitting.ChangeBeamDirection(dir, exNormal, loc, nextLoc);

		Point3f &inNormal = m_facets[track[i-1]].in_normal;
		p2 = ProjectPointToPlane(p1, dir, inNormal);
		double len = Length(p2 - p1);

		if (nextLoc == Location::In)
		{	// add internal path only
			m_splitting.ComputeCosA(dir, exNormal);
			double reRi = m_splitting.ComputeEffectiveReRi();
			len *= sqrt(reRi);
		}

#ifdef _DEBUG // DEB
		Point3f dddd = inNormal;
		dddd.d_param = -dddd.d_param;
		Point3f p22 = ProjectPointToPlane(p1, dir, dddd);
		double len1 = Length(p1 - p22);
		len1 *= sqrt(real(m_splitting.GetRi()));
		path1 += len1;
#endif
		path += len;

		p1 = p2;
		loc = nextLoc;
	}

#ifdef _DEBUG // DEB
//	path *= real(m_splitting.GetRi());
//	Point3f nFar1 = m_incidentDir;
//	Point3f nFar2 = -beam.direction;
//	double dd1 = m_splitting.FAR_ZONE_DISTANCE + DotProductD(p2, nFar1);
//	double dd2 = fabs(DotProductD(sourcePoint, nFar2) + m_splitting.FAR_ZONE_DISTANCE);
//	path += dd1;
//	path += dd2;
//	if (fabs(path - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}
