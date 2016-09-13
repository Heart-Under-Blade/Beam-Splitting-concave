#pragma once

#include <math.h>
#include "global.h"
#include "math/compl.hpp"

#define MAX_FACET_NUM 64

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{
public:
	Particle(double p_radius, double p_halfHeight, const complex &p_refractionIndex)
	{
		radius = p_radius;
		halfHeight = p_halfHeight;
		refractionIndex = p_refractionIndex;

		double re = real(refractionIndex);
		double im = imag(refractionIndex);
		double re_sqr =  re*re;
		refrI_sqr_re = re_sqr - im*im;
		refrI_sqr_im = 4*re_sqr*im;
	}

	virtual ~Particle() {}

	/**
	 * Rotate The angles of rotate are calculated by origin position
	 */
	virtual void Rotate(double beta, double gamma, double alpha) {}

public:
	Point3f facets[MAX_FACET_NUM][MAX_VERTEX_NUM];
	Point3f normals[MAX_FACET_NUM];					///< internal normals of facets
	Point3f externalNormals[MAX_FACET_NUM];			///< external normals of facets
	int facetNum;
	int vertexNums[MAX_FACET_NUM];

	complex refractionIndex;	///< complex value of refraction index of the particle
	double refrI_sqr_re;			///< real part of the sqr of refraction index
	double refrI_sqr_im;			///< imaginary part of the sqr of refraction index

protected:
	Point3f originNormals[MAX_FACET_NUM];
	Point3f originFacets[MAX_FACET_NUM][MAX_VERTEX_NUM];

	double radius;
	double halfHeight;
};

