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

Scattering::Scattering(const complex &refractiveIndex, int maxActNo,
					   double minBeamEnergy, double farFresnelZone)
	: m_maxActNo(maxActNo),
	  m_particle(nullptr),
	  m_regularIncidence(),
	  m_normalIncidence(),
	  m_totalReflectionIncidence(),
	  m_hasTracks(false),
	  m_minBeamEnergy(minBeamEnergy),
	  m_farFresnelZone(farFresnelZone)
{
	SetRefractiveIndex(refractiveIndex);

	// default start beam params
	Vector3f dir(0, 0, -1);
	Point3f point = dir * m_farFresnelZone;
	dir.d_param = Point3f::DotProduct(point, dir);
	Vector3f polarBasis(0, 1, 0);

	SetStartBeam(dir, polarBasis);
}

Scattering::~Scattering()
{
	delete m_shadowBeam.facet;
	delete m_splitting;
	delete m_startBeam.facet;
}

void Scattering::SetActiveParticle(Particle *particle)
{
	m_particle = particle;
	m_particle->GetPartByFacet(nullptr, m_workFacets);
}

void Scattering::SetStartBeam(const Vector3f &direction,
							  const Vector3f &polarBasis)
{
	m_startBeam.direction = direction;
	m_startBeam.direction.d_param = direction.d_param;
	m_startBeam.polarizationBasis = polarBasis;
	m_startBeam.facet = new Facet();
	m_startBeam.facet->index = -1;

	SetShadowBeamsParams();
}

void Scattering::SetShadowBeamsParams()
{
	Matrix2x2c jones;
	jones.m11 = -jones.m11;
	jones.m22 = -jones.m22;
	m_shadowBeam.Jones = jones;

	m_shadowBeam.direction = m_startBeam.direction;
	m_shadowBeam.polarizationBasis = m_startBeam.polarizationBasis;
	m_shadowBeam.facet = new Facet;
	m_shadowBeam.facet->index = INT_MAX;
}

void Scattering::Scatter(std::vector<Beam> *scatteredBeams)
{
	Reset();
	m_scatteredBeams = scatteredBeams;

	SplitStartBeam();
	SplitSecondaryBeams();
}

void Scattering::Scatter(TrackNode */*trackTree*/,
						 std::vector<Beam> &/*scatteredBeams*/)
{
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
		orderedPolygon.AddVertex(Point3f(base.x, base.y, m_farFresnelZone));
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

const Beam &Scattering::ExtractShadowBeam()
{
	std::vector<Point2f> projected;
	m_particle->ProjectToXYPlane(projected);

	std::vector<Point2f> projectedCleared;
	RemoveDublicatedVertices2f(projected, projectedCleared);

	OrderVertices2f(projectedCleared, m_shadowBeam);
	return m_shadowBeam;
}

void Scattering::SetIncidence()
{
	IncidenceType incType = m_splitting->GetIncidenceType();

	if (incType == IncidenceType::Regular)
	{
		m_incidence = &m_regularIncidence;
	}
	else if (incType == IncidenceType::Normal)
	{
		m_incidence = &m_normalIncidence;
	}
	else if (incType == IncidenceType::TotalReflection)
	{
		m_incidence = &m_totalReflectionIncidence;
	}
}

void Scattering::ComputeOpticalBeamParams(const Facet *facet, Beam beam,
										  const Polygon &resultShape)
{
	m_splitting->SetBeams(resultShape);
	auto &beams = m_splitting->beams;

	const Vector3f &normal = facet->normal[m_isBeamInside];
	m_splitting->ComputeSplittingParams(beam.direction, normal, m_isBeamInside);

	SetIncidence();

	m_incidence->ComputeDirections(beam, beams, m_isBeamInside);
	m_incidence->ComputeJonesMatrices(beam, beams, m_isBeamInside);
//	incidence->ComputeOpticalPaths(beam, beams);
}

Array<Facet*> *Scattering::FindVisibleFacets(const Beam &beam,
											 FacetChecker &checker,
											 Array<Facet*> &facets,
											 Array<Facet*> *visibleFacets)
{
	for (int i = 0; i < facets.nElems; ++i)
	{
		Facet *facet = facets.elems[i];

		if (checker.IsVisibleFacet(facet, beam, m_isBeamInside))
		{
			visibleFacets->Add(facet);
		}
	}
}

void Scattering::ComputeFacetEnergy(const Vector3f &facetNormal,
									const Polygon &lightedPolygon)
{
	double cosA = Point3f::DotProduct(m_startBeam.direction, facetNormal);
	m_incidentEnergy += lightedPolygon.Area() * cosA;
}

void Scattering::StackBeam(const Beam &beam)
{
#ifdef MODE_FIXED_OR
	beam.dirs.push_back(beam.direction);
	beam.pols.push_back(beam);
#endif
	m_internalBeams.Push(beam);
}

bool Scattering::IsFinalAct(const Beam &beam)
{
	return (beam.actNo >= m_maxActNo) ||
			(beam.Jones.Norm() < m_minBeamEnergy); // OPT: проверить сложность Norm()
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

const Beam &Scattering::GetStartBeam() const
{
	return m_startBeam;
}

Point3f Scattering::ComputeBeamDirection(const Vector3f &oldDir,
										 const Vector3f &normal,
										 bool isIn1, bool isIn2)
{
	Point3f newDir;

	if (!isIn1)
	{
		m_splitting->ComputeSplittingParams(oldDir, -normal, isIn1);
	}
	else
	{
		m_splitting->ComputeSplittingParams(oldDir, normal, isIn1);
	}

	if (isIn1 == isIn2)
	{
		m_splitting->ComputeReflectedDirection(newDir);
	}
	else
	{
		m_splitting->ComputeRefractedDirection(newDir);
	}

	return newDir;
}

OpticalPath Scattering::ComputeOpticalPath(const BeamInfo &info,
										   const Point3f &startPoint)
{
	OpticalPath path;

	Point3f dir = -info.beam->direction; // back direction
	bool isIn1 = false;
	bool isIn2;

	Point3f p1 = startPoint;
	Point3f p2;

	Facet *f1, *f2;

	auto &tr = info.track;

	// back tracing
	for (int act = tr.size()-1; act > 0; --act)
	{
		f1 = m_particle->GetActualFacet(tr[act]);
		f2 = m_particle->GetActualFacet(tr[act-1]);

		isIn2 = info.beam->IsInsideAtAct(act-1);
		dir = ComputeBeamDirection(dir, f1->ex_normal, isIn1, isIn2);
		p2 = Geometry::ProjectPointToPlane(p1, dir, f2->in_normal);
		double len = Point3f::Length(p2 - p1);

		if (isIn2)
		{
#ifdef _DEBUG // DEB
			len *= sqrt(real(m_splitting->GetRi()));
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
//	path *= real(m_splitting->GetRi());
//	Point3f nFar1 = m_originalBeam.direction;
//	Point3f nFar2 = -info.beam->direction;
//	double dd1 = m_splitting->FAR_ZONE_DISTANCE + Point3f::DotProduct(p2, nFar1);
//	double dd2 = fabs(Point3f::DotProduct(startPoint, nFar2) + m_splitting->FAR_ZONE_DISTANCE);

//	path.external += dd1;
//	path.external += dd2;

//	if (fabs(path.GetTotal() - beam.opticalPath) > 1)
//		int ff = 0;
#endif
	return path;
}

void Scattering::SelectOriginVisibleFacets(Array<Facet*> *facets)
{
	FindVisibleFacets(m_startBeam, m_lightChecker, m_workFacets, facets);
}

double Scattering::MeasureOpticalPath(const BeamInfo &info,
									  const Point3f sourcePoint)
{
}

double Scattering::MeasureFullOpticalPath(const BeamInfo &info,
										  const Point3f sourcePoint)
{
}

void Scattering::SetRefractiveIndex(const complex &ri)
{
	m_refractiveIndex = ri;
	m_splitting = new Splitting(ri);
	m_regularIncidence.SetSplitting(m_splitting);
	m_normalIncidence.SetSplitting(m_splitting);
	m_totalReflectionIncidence.SetSplitting(m_splitting);
}

void Scattering::SetFarFresnelZone(double value)
{
	m_farFresnelZone = value;
}

double Scattering::GetFarFresnelZone() const
{
	return m_farFresnelZone;
}

void Scattering::SplitBeamByVisibleFacets(Beam &beam)
{
	Array<Facet*> visibleFacets; // NOTE: DO NOT MAKE IT FIELD OF CLASS!!
	SelectVisibleFacets(beam, &visibleFacets);

	for (int i = 0; !isFinalFacet(i, visibleFacets); ++i)
	{
		Facet *facet = visibleFacets.elems[i];

		Polygon beamShape;
		bool isIntersected = Geometry::IncidentBeamToFacet(facet, beam, m_isBeamInside,
														   beam.direction, beamShape);
		if (isIntersected)
		{
			ComputeOpticalBeamParams(facet, beam, beamShape);
			ResolveBeams(beam, facet);
		}
	}
}

bool Scattering::isFinalFacet(int index, Array<Facet*> &facets)
{
	return index >= facets.nElems;
}

void Scattering::ResolveBeams(Beam &parentBeam, Facet *facet)
{
	auto &beams = m_splitting->beams;
	Track tr = parentBeam;
	tr.Update(facet);
	tr.RefreshId(facet->index, m_particle->nElems);

	beams.internal.CopyTrack(tr); // NOTE: this line must be before IsFinalAct

	if (!IsFinalAct(beams.internal))
	{
		beams.internal.SetIsInside(true);
		StackBeam(beams.internal);
	}

	if (m_splitting->type != IncidenceType::TotalReflection)
	{
		beams.external.CopyTrack(tr);
		beams.external.SetIsInside(false);
		ResolveExternalBeam(beams.external);
	}
}

void Scattering::Reset()
{
	m_isBeamInside = false;
	m_secondaryBeams = &m_internalBeams;
	m_internalBeams.Clear();
}
