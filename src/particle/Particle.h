#pragma once

#include <math.h>

#include "global.h"
#include "types.h"
#include "compl.hpp"
#include "vector_lib.h"

#define MAX_FACET_NUM 64

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

public:
	int facetNum;
	Point3f facets[MAX_FACET_NUM][MAX_VERTEX_NUM]; // REF: попробовать сделать структуру из facets и vertexNums
	int vertexNums[MAX_FACET_NUM];

	// REF: попробовать сделать структуру из 2 нормалей (снизит ли время трассировки?)
	Point3f normals[MAX_FACET_NUM];					///< internal normals of facets
	Point3f externalNormals[MAX_FACET_NUM];			///< external normals of facets

	// REF: сделать протектед
	complex refractionIndex;	///< complex value of refraction index of the particle
	double refrI_coef_re;		///< real part of the sqr of refraction index
	double refrI_coef_im;		///< imaginary part of the sqr of refraction index

	double halfHeight; // REF: сделать протектед

protected:
	Point3f m_originNormals[MAX_FACET_NUM];
	double m_rotMatrix[3][3];				///< rotation matrix for vertices

	double m_radius;

protected:
	virtual void SetOriginNormals() {}
	virtual void SetFacetParams() {}

	void Init(double p_radius, double p_halfHeight, const complex &p_refractionIndex);
	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void SetDParams();
	void SetExternalNormals();
};

