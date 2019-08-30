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
					   int maxActNo, const complex &refractiveIndex)
	: m_maxActNo(maxActNo),
	  m_particle(particle),
	  m_refractiveIndex(refractiveIndex),
	  m_splitting(refractiveIndex),
	  m_regularIncidence(),
	  m_normalIncidence(),
	  m_totalReflectionIncidence(),
	  m_hasTracks(false)
{
	m_particle->GetPartByFacet(nullptr, m_workFacets);
	CreateOriginalBeam(incidentLight.direction, incidentLight.polarizationBasis);

	Matrix2x2c jones;
	jones.m11 = -jones.m11;
	jones.m22 = -jones.m22;
	m_shadowBeam.SetMatrix(jones);

	m_shadowBeam.direction = incidentLight.direction;
	m_shadowBeam.polarizationBasis = incidentLight.polarizationBasis;
	m_shadowBeam.facet->index = INT_MAX;
}

Scattering::~Scattering()
{
	delete m_originalBeam.facet;
}

void Scattering::CreateOriginalBeam(const Vector3f &dir, const Vector3f &basis)
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
	m_treeSize = 0;
	m_scatteredBeams = &scatteredBeams;

	SplitOriginalBeam(scatteredBeams);
	SplitSecondaryBeams(scatteredBeams);
}

void Scattering::ScatterLight(TrackNode */*trackTree*/,
							  std::vector<Beam> &/*scatteredBeams*/)
{
}

void Scattering::SplitSecondaryBeams(std::vector<Beam> &scatteredBeams)
{
	while (m_treeSize != 0)
	{
		Beam beam = m_propagatingBeams[--m_treeSize];

		SplitBeamByVisibleFacets(beam);

		if (IsFinalAct(beam))
		{
			ReleaseBeam(beam);
		}
	}
}

void Scattering::OrderVertices2f(std::vector<Point2f> &vertices,
								 Polygon &orderedPolygon)
{
	orderedPolygon.Clear();

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
	n.d_param = m_particle->MaximalDimension();

	for (int i = 0; i < m_particle->nElems; i++)
	{
		auto &f = m_particle->elems[i].actual;

		if (Point3f::DotProduct(f.in_normal, -n) < EPS_COS_90)
		{
			for (int j = 0; j < f.nVertices; j++)
			{
//				auto p = ProjectPointToPlane(f.arr[j], -m_incidentDir, n);
//				projected.push_back(Point2f(-p.coordinates[0], -p.coordinates[1]));
				double tmp = (n.d_param - Point3f::DotProduct(n, f.vertices[j]));
				auto p = f.vertices[j] + n*tmp;
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

const Beam &Scattering::SetShadowBeam()
{
	std::vector<Point2f> projected;
	ProjectParticleToXY(projected);

	std::vector<Point2f> projectedCleared;
	RemoveDublicatedVertices2f(projected, projectedCleared);

	OrderVertices2f(projectedCleared, m_shadowBeam);
	return m_shadowBeam;
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
		m_incidence = &m_totalReflectionIncidence;
	}
}

bool Scattering::ComputeOpticalBeamParams(Facet *facet, Beam beam,
										  const Polygon &resultShape)
{
	m_splitting.SetBeams(resultShape);
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
								   Array<Facet*> &facets,
								   Array<Facet*> &visibleFacets)
{
	for (int i = 0; i < facets.nElems; ++i)
	{
		Facet *facet = facets.elems[i];

		if (checker.IsVisibleFacet(facet, beam))
		{
			visibleFacets.Add(facet);
		}
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
	if (Scattering::IsFinalAct(beam))
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
}

bool Scattering::IsFinalAct(const Beam &beam)
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

double Scattering::GetIncidentEnergy() const
{
	return m_incidentEnergy;
}

Point3f Scattering::ComputeBeamDirection(const Vector3f &oldDir,
										 const Vector3f &normal,
										 bool isIn1, bool isIn2)
{
	Point3f newDir;

	if (!isIn1)
	{
		m_splitting.ComputeSplittingParams(oldDir, -normal, isIn1);
	}
	else
	{
		m_splitting.ComputeSplittingParams(oldDir, normal, isIn1);
	}

	if (isIn1 == isIn2)
	{
		m_splitting.ComputeReflectedDirection(newDir);
	}
	else
	{
		m_splitting.ComputeRefractedDirection(newDir);
	}

	return newDir;
}

OpticalPath Scattering::ComputeOpticalPath(const Beam &beam,
										   const Point3f &startPoint,
										   std::vector<int> track)
{
	OpticalPath path;

	Point3f dir = -beam.direction; // back direction
	bool isIn1 = false;
	bool isIn2;

	Point3f p1 = startPoint;
	Point3f p2;

	Facet *f1, *f2;

	// back tracing
	for (int act = track.size()-1; act > 0; --act)
	{
		f1 = m_particle->GetActualFacet(track[act]);
		f2 = m_particle->GetActualFacet(track[act-1]);

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
	FindVisibleFacets(m_originalBeam, m_lightChecker, m_workFacets, facets);
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

	for (int i = 0; !isFinalFacet(i, visibleFacets); ++i)
	{
		Facet *facet = visibleFacets.elems[i];

		Polygon beamShape;
		bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, beam.isInside,
														   beam.direction, beamShape);
		if (isIntersected)
		{
			bool hasOutBeam = ComputeOpticalBeamParams(facet, beam, beamShape);
			PushBeamsToBuffer(beam, facet, hasOutBeam);
		}
	}
}

bool Scattering::isFinalFacet(int index, Array<Facet*> &facets)
{
	return index >= facets.nElems;
}

void Scattering::PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
								   bool hasOutBeam)
{
	auto &beams = m_splitting.beams;
	Track tr = parentBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(facet->index, m_particle->nElems);

	beams.internal.CopyTrack(tr);
	beams.internal.SetLocation(true);
	PushBeamToTree(beams.internal);

	if (hasOutBeam)
	{
		beams.external.CopyTrack(tr);
		beams.external.SetLocation(false);
		PushBeamToTree(beams.external);
	}
}
