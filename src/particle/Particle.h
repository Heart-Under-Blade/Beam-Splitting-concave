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
{	// REF: разработать нормальную классовую структуру
public:
	Particle();

	void Rotate(double beta, double gamma, double alpha);

	const double &GetMainSize() const;
	const double &GetSymmetryAngle() const;
	const complex &GetRefractionIndex() const;

	bool IsUnshadowedExternal(int facetId) const;
	bool IsShadowedInternal(int facetId) const;

public:
	Facet facets[MAX_FACET_NUM];	///< all facets of particle
	int m_facetNum;					///< number of facets

protected:
	struct ParticleState
	{
		Facet facets[MAX_FACET_NUM];
	};

	ParticleState defaultState;

	double m_mainSize;			///< max size of particle (diameter or height or smth)
	double m_symmetryAngle;		///< angle of particle symmetry
	complex m_refractionIndex;	///< complex value of refraction index of the particle

	IntArray m_unshadowedExternalFacets;
	IntArray m_shadowedInternalFacets;

protected:
	void Init(int m_facetNum, const complex &refrIndex, double symAngle,
			  double size);

	void SetDefaultNormals();
	void SetActualState();
	virtual void SetFacetParams() {}

	void SetRotateMatrix(double beta, double gamma, double alpha);
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);

private:
	void SetDParams();

private:
	double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
};

