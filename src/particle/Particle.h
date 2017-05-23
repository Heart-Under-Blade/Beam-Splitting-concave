#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"

#define ROT_MTR_RANK 3

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle
{
public:
	Particle();

	void Rotate(double beta, double gamma, double alpha);

	const double &GetMainSize() const;
	const double &GetSymmetryBeta() const;
	const double &GetSymmetryGamma() const;
	const complex &GetRefractionIndex() const;

public:
	Facet facets[MAX_FACET_NUM];	///< all facets of particle
	int facetNum;					///< number of facets

protected:
	Facet defaultFacets[MAX_FACET_NUM];

	double m_mainSize;			///< max size of particle (diameter or height or smth)
	double m_symmetryBeta;		///< angle of particle symmetry of beta angle
	double m_symmetryGamma;		///< angle of particle symmetry of gamma angle
	complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
	void Init(int facetCount, const complex &refrIndex,
			  double symGamma, double symBeta, double size);

	void SetDefaultNormals();
	void SetDefaultCenters();
	void SetActualState();
	virtual void SetFacetParams() {}

	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);

private:
	void SetDParams();

private:
	double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
	void RotateCenters();
};

