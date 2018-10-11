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
	  m_nActsMax(nActs)
{
	m_originBeam.direction = incidentLight->direction;
	m_originBeam.direction.d_param = incidentLight->direction.d_param;
	m_originBeam.polarizationBasis = incidentLight->polarizationBasis;
	m_originBeam.isInside = false;

	m_splitting.ComputeRiParams(m_particle->GetRefractiveIndex());
}

void Scattering::ScatterLight(std::vector<Beam> &scatteredBeams)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
#endif
	m_treeSize = 0;

	SplitLightToBeams(scatteredBeams);
	SplitBeams(scatteredBeams);
}

// REF: перенести в Tracks
IdType Scattering::RecomputeTrackId(const IdType &oldId, int facetId)
{
	return (oldId + (facetId + 1)) * (m_particle->nElems + 1);
}

void Scattering::PushBeamToTree(Beam &beam, Facet *facet, int level, bool isIn)
{
	beam.SetTracingParams(facet, level, isIn);
#ifdef _DEBUG // DEB
	beam.dirs.push_back(beam.direction);
#endif
	m_tracingBeams[m_treeSize++] = beam;
}

bool Scattering::ComputeOpticalBeamParams(Facet *facet, Beam beam)
{
#ifdef _DEBUG // DEB
	m_splitting.inBeam.pols = beam.pols;
	m_splitting.outBeam.pols = beam.pols;
#endif
	m_splitting.ComputeParams(beam.direction, facet->normal[!beam.isInside], beam.isInside);

	Incidence *incidence = m_splitting.GetIncidence();
	incidence->ComputeDirections(beam, m_splitting);
	incidence->ComputeJonesMatrices(beam, m_splitting);
	incidence->ComputeOpticalPaths(beam, m_splitting);

	return m_splitting.HasOutBeam();
}

Particle *Scattering::GetParticle() const
{
	return m_particle;
}

void Scattering::FindVisibleFacets(const Beam &beam, FacetChecker &checker,
								   int begin, int end, Array<Facet*> &facets)
{
	for (int i = begin; i < end; ++i)
	{
		Facet *facet = m_particle->GetActualFacet(i);

		if (checker.IsVisibleFacet(facet, beam))
		{
			facets.Add(facet);
		}
	}
}

void Scattering::ComputeFacetEnergy(const Vector3f &facetNormal,
									const Polygon &lightedPolygon)
{
	double cosA = Point3f::DotProduct(m_originBeam.direction, facetNormal);
	m_incidentEnergy += lightedPolygon.Area() * cosA;
}

void Scattering::PushBeamToTree(Beam &beam, const Beam &oldBeam,
								const IdType &newId, Facet *facet, bool isIn)
{
	beam.id = newId;
	beam.locations = oldBeam.locations;
#ifdef _DEBUG // DEB
	beam.dirs = oldBeam.dirs;
#endif
	PushBeamToTree(beam, facet, oldBeam.act+1, isIn);
}

void Scattering::RotateParticle(const Angle &angle)
{
	m_particle->Rotate(angle);
}

bool Scattering::IsTerminalAct(const Beam &beam)
{
	return (beam.act >= m_nActsMax) || (beam.J.Norm() < EPS_BEAM_ENERGY);
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

Point3f Scattering::ComputeBeamDirection(const Vector3f &oldDir,
										 const Vector3f &normal,
										 bool isIn1, bool isIn2)
{
	Point3f newDir;

	if (!isIn1)
	{
		m_splitting.ComputeParams(oldDir, -normal, isIn1);
	}
	else
	{
		m_splitting.ComputeParams(oldDir, normal, isIn1);
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
	Point3f nFar1 = m_originBeam.direction;
	Point3f nFar2 = -beam.direction;
	double dd1 = m_splitting.FAR_ZONE_DISTANCE + Point3f::DotProduct(p2, nFar1);
	double dd2 = fabs(Point3f::DotProduct(startPoint, nFar2) + m_splitting.FAR_ZONE_DISTANCE);

	path.external += dd1;
	path.external += dd2;

//	if (fabs(path.GetTotal() - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}
