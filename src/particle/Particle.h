#pragma once

#include <math.h>

#include "types.h"
#include "compl.hpp"
#include "geometry_lib.h"

#define in_normal normal[0]
#define ex_normal normal[1]

enum class Location: bool
{
	Internal, External
};

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{	// TODO: разработать нормальную классовую структуру
public:
	Particle();
	Particle(double radius, double halfHeight, const complex &refractionIndex);
	virtual ~Particle() {}

	/**
	 * Rotate The angles of rotate are calculated by origin position
	 */
	virtual void Rotate(double /*beta*/, double /*gamma*/, double /*alpha*/) {}

	double GetHalfHeight() const;
	const complex &GetRefractionIndex() const;

public:
	Facet facets[MAX_FACET_NUM];
	int facetNum;
	Point3f centers[MAX_FACET_NUM];

protected:
	Point3f m_originCenters[MAX_FACET_NUM];
	Point3f m_originNormals[MAX_FACET_NUM];

	double m_rotMatrix[3][3];	///< rotation matrix for vertices

	double m_radius;
	double m_halfHeight;

	complex m_refractionIndex;	///< complex value of refraction index of the particle

protected:
	virtual void SetOriginNormals() {}
	virtual void SetFacetParams() {}

	void Init(double radius, double halfHeight, const complex &refractionIndex);
	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void SetDParams();
	void SetExternalNormals();
	void CopyFacet(Point3f *points, Facet &result);
};

