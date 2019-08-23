#include "Scattering.h"

#include <float.h>
#include <assert.h>
#include <algorithm>

#include "macro.h"
#include "geometry_lib.h"

#ifdef _DEBUG // DEB
//#include <iostream>
#endif

#define NORM_CEIL	FLT_EPSILON + 1

using namespace std;

Scattering::Scattering(Particle *particle, const Light &incidentLight,
					   int maxActNo)
	: m_particle(particle),
	  m_maxActNo(maxActNo),
	  m_splitting(m_particle->GetRefractiveIndex()),
	  m_completeReflectionIncidence(),
	  m_regularIncidence(),
	  m_normalIncidence()
{
	CreateOriginBeam(incidentLight.direction, incidentLight.polarizationBasis);
}

Scattering::~Scattering()
{
	delete m_originalBeam.facet;
}

void Scattering::CreateOriginBeam(const Vector3f &dir, const Vector3f &basis)
{
	m_originalBeam.direction = dir;
	m_originalBeam.direction.d_param = dir.d_param;
	m_originalBeam.polarizationBasis = basis;
	m_originalBeam.isInside = false;
	m_originalBeam.facet = new Facet();
	m_originalBeam.facet->index = -1;
}

void Scattering::ScatterLight(std::vector<Beam> &scatteredBeams)
{
<<<<<<< HEAD
	beam.SetTracingParams(facetId, level, location);
#ifdef _DEBUG // DEB
//	beam.dirs.push_back(beam.direction);
=======
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#endif
	m_treeSize = 0;
	m_scatteredBeams = &scatteredBeams;

	SplitOriginalBeam(scatteredBeams);
	SplitSecondaryBeams(scatteredBeams);
}

void Scattering::SplitSecondaryBeams(std::vector<Beam> &scatteredBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_propagatingBeams[--m_treeSize];

		SplitBeamByVisibleFacets(beam);

		if (IsTerminalAct(beam))
		{
			ReleaseBeam(beam);
		}
	}
}

void Scattering::SetIncidence()
{
	IncidenceType incType = m_splitting.GetIncidenceType();

	if (incType == IncidenceType::Regular)
	{
		m_incidence = &m_regularIncidence;
	}
	else if (incType == IncidenceType::Normal)
	{
		m_incidence = &m_normalIncidence;
	}
	else if (incType == IncidenceType::CompleteReflection)
	{
		m_incidence = &m_completeReflectionIncidence;
	}
}

bool Scattering::ComputeOpticalBeamParams(Facet *facet, Beam beam)
{
	auto &beams = m_splitting.beams;

	const Vector3f &normal = facet->normal[beam.isInside];
	m_splitting.ComputeSplittingParams(beam.direction, normal, beam.isInside);

	SetIncidence();

	m_incidence->ComputeDirections(beam, beams);
	m_incidence->ComputeJonesMatrices(beam, beams);
//	incidence->ComputeOpticalPaths(beam, beams);

	return m_splitting.HasOutBeam();
}

void Scattering::FindVisibleFacets(const Beam &beam, FacetChecker &checker,
								   int begin, int end, Array<Facet*> &facets)
{
	for (int i = begin; i < end; ++i)
	{
<<<<<<< HEAD
		Point3f p = inBeam.Center();
		double path = m_splitting.ComputeIncidentOpticalPath(m_incidentDir, p);
		inBeam.opticalPath = 0;
		outBeam.opticalPath = 0;
		inBeam.AddOpticalPath(path);
		outBeam.AddOpticalPath(path);
		outBeam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(outBeam);
=======
		Facet *facet = m_particle->GetActualFacet(i);

		if (checker.IsVisibleFacet(facet, beam))
		{
			facets.Add(facet);
		}
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	}
}

void Scattering::ComputeFacetEnergy(const Vector3f &facetNormal,
									const Polygon &lightedPolygon)
{
	double cosA = Point3f::DotProduct(m_originalBeam.direction, facetNormal);
	m_incidentEnergy += lightedPolygon.Area() * cosA;
}

void Scattering::PushBeamToTree(Beam &beam)
{
<<<<<<< HEAD
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
=======
	if (Scattering::IsTerminalAct(beam))
	{
		if (!beam.isInside)
		{
			ReleaseBeam(beam);
		}
	}
	else
	{
#ifdef MODE_FIXED_OR
		beam.dirs.push_back(beam.direction);
		beam.pols.push_back(beam);
#endif
		m_propagatingBeams[m_treeSize++] = beam;
	}
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
	return (beam.actNo >= m_maxActNo) || (beam.Jones.Norm() < EPS_BEAM_ENERGY);
}

/** NOTE: Result beams are ordered in inverse direction */
void Scattering::SetPolygonByFacet(Facet *facet, Polygon &polygon) const
{
	int size = facet->nVertices;
	polygon.nVertices = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.vertices[i] = facet->vertices[size-i];
	}
}

<<<<<<< HEAD
	__m128 _clip_normal = _mm_setr_ps(clipNormal.cx, clipNormal.cy, clipNormal.cz, 0.0);

	int clipSize = clip.nVertices;
	__m128 _diff_pol[MAX_VERTEX_NUM];
=======
double Scattering::GetIncidentEnergy() const
{
	return m_incidentEnergy;
}
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

Point3f Scattering::ComputeBeamDirection(const Vector3f &oldDir,
										 const Vector3f &normal,
										 bool isIn1, bool isIn2)
{
	Point3f newDir;

<<<<<<< HEAD
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
=======
	if (!isIn1)
	{
		m_splitting.ComputeSplittingParams(oldDir, -normal, isIn1);
	}
	else
	{
		m_splitting.ComputeSplittingParams(oldDir, normal, isIn1);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	}

	if (isIn1 == isIn2)
	{
		m_splitting.ComputeReflectedDirection(newDir);
	}
<<<<<<< HEAD

	for (int i = 0; i < polygon.nVertices; ++i)
=======
	else
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	{
		m_splitting.ComputeRefractedDirection(newDir);
	}

	return newDir;
}

OpticalPath Scattering::ComputeOpticalPath(const Beam &beam,
										   const Point3f &startPoint,
										   std::vector<int> track)
{
<<<<<<< HEAD
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
=======
	OpticalPath path;

	Point3f dir = -beam.direction; // back direction
	bool isIn1 = false;
	bool isIn2;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	Point3f p1 = startPoint;
	Point3f p2;

<<<<<<< HEAD
	Point3f p2 = m_particle->facets[facetID].arr[facetSize-1];
	_p2 = _mm_load_ps(p2.coordinates);
=======
	Facet *f1, *f2;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	// back tracing
	for (int act = track.size()-1; act > 0; --act)
	{
<<<<<<< HEAD
		_p1 = _p2;
		p2 = m_particle->facets[facetID].arr[i];
		_p2 = _mm_load_ps(p2.coordinates);

		bufferSize = outputSize;
		outputSize = 0;
=======
		f1 = m_particle->GetActualFacet(track[act]);
		f2 = m_particle->GetActualFacet(track[act-1]);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

		isIn2 = beam.IsInsideOnAct(act-1);
		dir = ComputeBeamDirection(dir, f1->ex_normal, isIn1, isIn2);
		p2 = Geometry::ProjectPointToPlane(p1, dir, f2->in_normal);
		double len = Point3f::Length(p2 - p1);

		if (isIn2)
		{
#ifdef _DEBUG // DEB
			len *= sqrt(real(m_splitting.GetRi()));
#endif
			path.internal += len;
		}
		else
		{
			path.external += len;
		}

		p1 = p2;
		isIn1 = isIn2;
	}

#ifdef _DEBUG // DEB
//	path *= real(m_splitting.GetRi());
	Point3f nFar1 = m_originalBeam.direction;
	Point3f nFar2 = -beam.direction;
//	double dd1 = m_splitting.FAR_ZONE_DISTANCE + Point3f::DotProduct(p2, nFar1);
//	double dd2 = fabs(Point3f::DotProduct(startPoint, nFar2) + m_splitting.FAR_ZONE_DISTANCE);

//	path.external += dd1;
//	path.external += dd2;

//	if (fabs(path.GetTotal() - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}

void Scattering::SelectOriginVisibleFacets(Array<Facet*> &facets)
{
	FindVisibleFacets(m_originalBeam, m_lightChecker, 0, m_particle->nElems, facets);
}

void Scattering::ReleaseBeam(Beam &beam)
{
//	beam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
	m_scatteredBeams->push_back(beam);
}

void Scattering::SplitBeamByVisibleFacets(Beam &beam)
{
	Array<Facet*> visibleFacets; // NOTE: DO NOT MAKE IT FIELD OF CLASS!!
	SelectVisibleFacets(beam, visibleFacets);

	for (int i = 0; !isTerminalFacet(i, visibleFacets); ++i)
	{
		Facet *facet = visibleFacets.elems[i];

		Polygon beamShape;
		bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, beam.isInside,
														   beam.direction, beamShape);
		if (isIntersected)
		{
<<<<<<< HEAD
			p.cx = _output_points[i][0];
			p.cy = _output_points[i][1];
			p.cz = _output_points[i][2];
			polygon.arr[polygon.nVertices++] = p;
=======
			m_splitting.SetBeams(beamShape);
			bool hasOutBeam = ComputeOpticalBeamParams(facet, beam);
			PushBeamsToBuffer(beam, facet, hasOutBeam);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
		}
	}
}

bool Scattering::isTerminalFacet(int index, Array<Facet*> &facets)
{
<<<<<<< HEAD
	const Polygon &facet = m_facets[facetId];
	int size = facet.nVertices;
	polygon.nVertices = size;
	--size;

	for (int i = 0; i <= size; ++i)
	{
		polygon.arr[i] = facet.arr[size-i];
	}
=======
	return index >= facets.nElems;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}

void Scattering::PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
								   bool hasOutBeam)
{
<<<<<<< HEAD
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
//		Point3f dddd = inNormal;
//		dddd.d_param = -dddd.d_param;
//		Point3f p22 = ProjectPointToPlane(p1, dir, dddd);
//		double len1 = Length(p1 - p22);
//		len1 *= sqrt(real(m_splitting.GetRi()));
//		path1 += len1;
=======
	auto &beams = m_splitting.beams;
	Track tr = parentBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(facet->index, m_particle->nElems);

#ifdef _DEBUG // DEB
	if (tr.id == 4222009)
		int fff = 0;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#endif
	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);
	PushBeamToTree(beams.internal);

	if (hasOutBeam)
	{
		beams.external.CopyTrack(tr);
		beams.external.SetLocation(false);
		PushBeamToTree(beams.external);
	}
<<<<<<< HEAD

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
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}
