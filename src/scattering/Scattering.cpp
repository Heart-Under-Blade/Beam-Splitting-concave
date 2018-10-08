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
	  m_maxNActs(nActs)
{
	m_originBeam.direction = incidentLight->direction;
	m_originBeam.direction.d_param = incidentLight->direction.d_param;
	m_originBeam.polarizationBasis = incidentLight->polarizationBasis;
	m_originBeam.isInside = false;

	m_splitting.ComputeRiParams(m_particle->GetRefractiveIndex());
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

bool Scattering::SetOpticalBeamParams(Facet *facet, Beam beam)
{
#ifdef _DEBUG // DEB
	m_splitting.inBeam.pols = beam.pols;
	m_splitting.outBeam.pols = beam.pols;
#endif
	bool hasOutBeam = true;

	if (m_splitting.IsNormalIncidence()) // normal incidence
	{
		ComputeOpticalParams(m_normalIncidence, beam);
	}
	else // regular incidence
	{
		if (beam.isInside)
		{
			m_splitting.ComputeSplittingParams(beam.direction, facet->ex_normal, beam.isInside);
			ComputePolarisationParams(beam.direction, facet->ex_normal, beam);

			hasOutBeam = !m_splitting.IsCompleteReflection();

			if (hasOutBeam)
			{
				ComputeOpticalParams(m_regularIncidence, beam);
			}
			else // complete internal reflection incidence
			{
				ComputeOpticalParams(m_completeReflectionIncidence, beam);
			}
		}
		else // beam is external
		{
			m_splitting.ComputeCosA(facet->in_normal, beam.direction);
			m_splitting.ComputeSplittingParams(beam.direction, facet->in_normal, beam.isInside);
			ComputePolarisationParams(-beam.direction, facet->in_normal, beam);
			ComputeOpticalParams(m_regularIncidence, beam);
		}
	}

	return hasOutBeam;
}

void Scattering::ComputeOpticalParams(const Incidence &incidence,
									  const Beam &incidentBeam)
{
	incidence.ComputeDirections(incidentBeam, m_splitting);
	incidence.ComputeJonesMatrices(incidentBeam, m_splitting);
	incidence.ComputeOpticalPaths(incidentBeam, m_splitting);
}

void Scattering::ComputePolarisationParams(const Vector3f &dir,
										   const Vector3f &facetNormal, Beam &beam)
{
	Point3f newBasis = Point3f::CrossProduct(facetNormal, dir);
	Point3f::Normalize(newBasis);
	beam.RotateJMatrix(newBasis);
	beam.polarizationBasis = newBasis;
}

void Scattering::SplitLightToBeams(Facet *facet)
{
}

Particle *Scattering::GetParticle() const
{
	return m_particle;
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

// TODO: пофиксить
void Scattering::ScatterLight(const std::vector<std::vector<int>> &/*tracks*/,
							  std::vector<Beam> &/*outBeams*/)
{
//	m_particle->Rotate(beta, gamma, 0);

//	for (int i = 0; i < tracks.size(); ++i)
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

//		int size = tracks.at(i).size();

//		try // internal beams
//		{
//			for (int j = 1; j < size; ++j)
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
	return (beam.act >= m_maxNActs) || (beam.J.Norm() < EPS_BEAM_ENERGY);
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
		m_splitting.ComputeCosA(oldDir, -normal);
		m_splitting.ComputeSplittingParams(oldDir, -normal, isIn1);
	}
	else
	{
		m_splitting.ComputeCosA(oldDir, normal);
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
			m_splitting.ComputeCosA(dir, f1->ex_normal);
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
