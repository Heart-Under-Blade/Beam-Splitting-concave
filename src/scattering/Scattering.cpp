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

Scattering::Scattering(Particle *particle, const Light &incidentLight,
					   int maxActNo, const complex &refractiveIndex)
	: m_particle(particle),
	  m_maxActNo(maxActNo),
	  m_splitting(refractiveIndex),
	  m_completeReflectionIncidence(),
	  m_regularIncidence(),
	  m_refractiveIndex(refractiveIndex),
	  m_normalIncidence(),
	  m_hasTracks(false)
{
	m_particle->GetPartByFacet(nullptr, m_workFacets);
	CreateOriginalBeam(incidentLight.direction, incidentLight.polarizationBasis);
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

		if (m_hasTracks)
		{
			m_trackTreeNode = m_trackTreeNode->FindNode(beam.facet->index);
		}

		SplitBeamByVisibleFacets(beam);

		if (IsFinalAct(beam))
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

#ifdef _DEBUG // DEB
	if (tr.id == 4222009)
		int fff = 0;
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
}
