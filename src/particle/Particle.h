#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"

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

	// TODO: допилить
	void SetFromFile(const char *filename);

	void Rotate(double beta, double gamma, double alpha);

	const double &GetMainSize() const;
	const complex &GetRefractionIndex() const;
	const Symmetry &GetSymmetry() const;
	virtual void GetAggPartFacetIDRange(int /*id*/, int &/*begin*/, int &/*end*/) const {}

	virtual bool IsComlicated() const;

	void Output();

public:
	Facet facets[MAX_FACET_NUM];	///< all facets of particle
	int facetNum;					///< number of facets
	bool isAggregate = false;

protected:
	Facet defaultFacets[MAX_FACET_NUM];

	double m_mainSize;			///< max size of particle (diameter or height or smth)	
	Symmetry m_symmetry;		///< angle of particle symmetry

	complex m_refractiveIndex;	///< complex value of refractive index of the particle

protected:
	void Init(int facetCount, const complex &refrIndex, double size);

	void SetDefaultNormals();
	void SetDefaultCenters();
	void SetActualState();
	void SetSymmetry(double beta, double gamma, double alpha = 0);
	virtual void SetFacetParams() {}

private:
	void SetDParams();
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void RotateCenters();
	void SetRotateMatrix(double beta, double gamma, double alpha);

private:
	double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
};

