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
					   int maxActNo)
	: m_particle(particle),
	  m_maxActNo(maxActNo),
	  m_splitting(m_particle->GetRefractiveIndex())
{
	CreateOriginBeam(incidentLight.direction, incidentLight.polarizationBasis);
}

Scattering::~Scattering()
{
	delete m_originBeam.facet;
}

void Scattering::CreateOriginBeam(const Vector3f &dir, const Vector3f &basis)
{
	m_originBeam.direction = dir;
	m_originBeam.direction.d_param = dir.d_param;
	m_originBeam.polarizationBasis = basis;
	m_originBeam.isInside = false;
	m_originBeam.facet = new Facet();
	m_originBeam.facet->index = -1;
}

void Scattering::ScatterLight(std::vector<Beam> &scatteredBeams)
{
#ifdef _CHECK_ENERGY_BALANCE
	m_incidentEnergy = 0;
#endif
	m_treeSize = 0;
	m_scatteredBeams = &scatteredBeams;

	SplitOriginBeam(scatteredBeams);
	SplitSecondaryBeams(scatteredBeams);
}

void Scattering::SplitSecondaryBeams(std::vector<Beam> &scatteredBeams)
{
#ifdef _DEBUG // DEB
	int count = 0;
#endif
	while (m_treeSize != 0)
	{
		Beam beam = m_propagatingBeams[--m_treeSize];
#ifdef _DEBUG // DEB
		++count;
		if (beam.id == 1142698)
			int ffgf = 0;
//		cout << count << endl;
#endif
		SplitBeamByVisibleFacets(beam);

		if (IsTerminalAct(beam))
		{
			ReleaseBeam(beam);
		}
	}
}

bool Scattering::ComputeOpticalBeamParams(Facet *facet, Beam beam)
{
#ifdef _DEBUG // DEB
	m_splitting.inBeam.pols = beam.pols;
	m_splitting.outBeam.pols = beam.pols;
#endif
	const Vector3f &normal = facet->normal[beam.isInside];
	m_splitting.ComputeParams(beam.direction, normal, beam.isInside);

	Incidence *incidence = m_splitting.GetIncidence();
	incidence->ComputeDirections(beam, m_splitting);
	incidence->ComputeJonesMatrices(beam, m_splitting);
	incidence->ComputeOpticalPaths(beam, m_splitting);

	return m_splitting.HasOutBeam();
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

void Scattering::PushBeamToTree(Beam &beam)
{
	if (Scattering::IsTerminalAct(beam))
	{
		if (!beam.isInside)
		{
			ReleaseBeam(beam);
		}
	}
	else
	{
#ifdef _DEBUG // DEB
		beam.dirs.push_back(beam.direction);
#endif
		m_propagatingBeams[m_treeSize++] = beam;
	}
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
		polygon.arr[i] = facet->arr[size-i];
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

void Scattering::ReleaseBeam(Beam &beam)
{
	beam.opticalPath += m_splitting.ComputeOutgoingOpticalPath(beam); // добираем оптический путь
	m_scatteredBeams->push_back(beam);
}

void Scattering::SplitBeamByVisibleFacets(Beam &beam)
{
	Array<Facet*> visibleFacets; // NOTE: DO NOT MAKE IT FIELD OF CLASS!!
	SelectVisibleFacets(beam, visibleFacets);

	for (int i = 0; !isTerminalFacet(i, visibleFacets); ++i)
	{
		Facet *facet = visibleFacets.elems[i];
#ifdef _DEBUG // DEB
		if (i == 9)
			int fff = 0;
#endif
		Polygon beamShape;
		bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, beam.isInside,
														   beam.direction, beamShape);
		if (isIntersected)
		{
			m_splitting.SetBeams(beamShape);
			bool hasOutBeam = ComputeOpticalBeamParams(facet, beam);
			PushBeamsToBuffer(beam, facet, hasOutBeam);
		}
	}
}

bool Scattering::isTerminalFacet(int index, Array<Facet*> &facets)
{
	return index >= facets.nElems;
}

void Scattering::PushBeamsToBuffer(Beam &parentBeam, Facet *facet,
								   bool hasOutBeam)
{
	Track tr = parentBeam;
	tr.Update(facet);
	tr.RecomputeTrackId(facet->index, m_particle->nElems);

#ifdef _DEBUG // DEB
	if (tr.id == 4222009)
		int fff = 0;
#endif
	m_splitting.inBeam.CopyTrack(tr);
	m_splitting.inBeam.SetLocation(true);
#ifdef _DEBUG // DEB
	m_splitting.inBeam.dirs = parentBeam.dirs;
#endif
	PushBeamToTree(m_splitting.inBeam);

	if (hasOutBeam)
	{
		m_splitting.outBeam.CopyTrack(tr);
		m_splitting.outBeam.SetLocation(false);
#ifdef _DEBUG // DEB
		m_splitting.outBeam.dirs = parentBeam.dirs;
#endif
		PushBeamToTree(m_splitting.outBeam);
	}
}
