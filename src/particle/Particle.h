#pragma once

#include <math.h>

#include "types.h"
#include "compl.hpp"
#include "vector_lib.h"

#define MAX_FACET_NUM 64

#define MIN_VERTEX_NUM 3
#define MAX_VERTEX_NUM 64

#define in_normal normal[0]
#define ex_normal normal[1]

enum class Location: bool
{
	Internal, External
};

struct Polygon
{
	Point3f arr[MAX_VERTEX_NUM];
	int size = 0;
};

struct Facet
{
	Polygon polygon;
	Point3f normal[2]; ///< internal and external
};

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{
public:
	Particle();
	Particle(double p_radius, double p_halfHeight, const complex &p_refractionIndex);
	virtual ~Particle() {}

	/**
	 * Rotate The angles of rotate are calculated by origin position
	 */
	virtual void Rotate(double /*beta*/, double /*gamma*/, double /*alpha*/) {}

	double GetHalfHeight() const;

public:
	Facet facets[MAX_FACET_NUM];
	int facetNum;

	// REF: сделать протектед
	complex refractionIndex;	///< complex value of refraction index of the particle
	double ri_coef_re;			///< real part of the sqr of refraction index
	double ri_coef_im;			///< imaginary part of the sqr of refraction index

	Point3f centers[MAX_FACET_NUM];

protected:
	Point3f m_originCenters[MAX_FACET_NUM];
	Point3f m_originNormals[MAX_FACET_NUM];

	double m_rotMatrix[3][3];				///< rotation matrix for vertices

	double m_radius;
	double m_halfHeight;

protected:
	virtual void SetOriginNormals() {}
	virtual void SetFacetParams() {}

	void Init(double p_radius, double p_halfHeight, const complex &p_refractionIndex);
	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void SetDParams();
	void SetExternalNormals();
	void CopyFacet(Point3f *points, Facet &result);
};

