#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{	// REF: разработать нормальную классовую структуру
public:
	Particle();
	Particle(double radius, double halfHeight, const complex &refractionIndex);

	virtual void Rotate(double beta, double gamma, double alpha);

	const double &GetHalfHeight() const;
	const complex &GetRefractionIndex() const;

	bool IsUnshadowedExternal(int facetId) const;
	bool IsShadowedInternal(int facetId) const;

public:
	Facet facets[MAX_FACET_NUM];	///< all facets of particle
	int facetNum;					///< number of facets
	Point3f centers[MAX_FACET_NUM]; ///< centers of facets (for fast access without calc)

protected:
	Point3f m_originCenters[MAX_FACET_NUM];
	Point3f m_originNormals[MAX_FACET_NUM];

	double m_rotMatrix[3][3];	///< rotation matrix for vertices

	double m_radius;
	double m_halfHeight;

	complex m_refractionIndex;	///< complex value of refraction index of the particle

	IntArray m_unshadowedExternalFacets;
	IntArray m_shadowedInternalFacets;

protected:
	virtual void SetDefaultNormals() {}
	virtual void SetFacetParams() {}

	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void CopyFacet(Point3f *points, Facet &result);
	void SetActualNormals();

private:
	void Init(double radius, double halfHeight, const complex &refractionIndex);
	void SetDParams();
	void SetExternalNormals();
	void SetInternalNormals();
};

