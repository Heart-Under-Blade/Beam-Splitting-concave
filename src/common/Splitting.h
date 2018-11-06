#pragma once

#include "Beam.h"

#include "RegularIncidence.h"
#include "NormalIncidence.h"
#include "CompleteReflectionIncidence.h"

#define EPS_ORTO_FACET 0.0001

template <class T>
class SplittedBeams
{
public:
	T internal;
	T external;
};

class Splitting
{
public:
	Splitting(const complex &ri);

	Incidence *GetIncidence() const;

	void ComputeParams(const Vector3f &dir, const Vector3f &normal,
					   bool isInside);

	void ComputePolarisationParams(Beam &beam);

	void ComputeReflectedDirection(Vector3f &dir) const
	{
		dir = r - m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

	void ComputeRefractedDirection(Vector3f &dir) const
	{
		dir = r/sqrt(s) + m_normal;
		Point3f::Normalize(dir); // REF, OPT: нужно ли это нормализовать всегда?
	}

	void SetBeams(const Polygon &beamShape);
	void SetNormal(const Point3f &normal);

	bool HasOutBeam();
	double ComputeEffectiveReRi() const;

	complex GetRi() const;

public:
	double cosA;
	Point3f r;
	double reRiEff;
	double s;

	SplittedBeams<Beam> beams;

	complex m_ri;	///< Refractive index of a Particle

	Point3f m_normal;
	Incidence *m_incidence;

private:
	bool m_hasOutBeam;

	double m_cRiRe;
	double m_cRiRe2;
	double m_cRiIm;

	RegularIncidence				m_regularIncidence;
	NormalIncidence					m_normalIncidence;
	CompleteReflectionIncidence		m_completeReflectionIncidence;
};


class FacetChecker
{
public:
	virtual bool IsVisibleFacet(Facet *facet, const Beam &beam)
	{
		if (facet->index != beam.facet->index)
		{
			const Point3f &facetNormal = facet->normal[beam.isInside];
			double cosA = Point3f::DotProduct(facetNormal, beam.direction);
			return cosA > FLT_EPSILON;
		}
		else
		{
			return false;
		}
	}
};

class LightFacetChecker : public FacetChecker
{
public:
};

class BeamFacetChecker : public FacetChecker
{
public:
	bool IsVisibleFacet(Facet *facet, const Beam &beam) override
	{
		if (FacetChecker::IsVisibleFacet(facet, beam))
		{
			Point3f vectorFromBeamToFacet = facet->center - beam.facet->center;
			const Point3f &beamNormal = beam.facet->normal[!beam.isInside];
			double cosBF = Point3f::DotProduct(beamNormal, vectorFromBeamToFacet);
			return (cosBF >= EPS_ORTO_FACET);
		}
		else
		{
			return false;
		}
	}
};
