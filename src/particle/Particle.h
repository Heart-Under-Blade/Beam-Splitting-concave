#pragma once

#include <math.h>
#include <vector>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"
#include "Rotator.h"

struct ParticleFacet
{
	Facet origin; // facet with origin coordinates of points
	Facet actual; // facet with actual coordinates of points
};

/**
 * @brief The Particle class is the base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle : public Array<ParticleFacet>
{
public:
	Particle();
	Particle(int nFacets, const complex &refrIndex, bool isNonConvex = false);

	Facet *GetActualFacet(int i);
	void SetFromFile(const std::string &filename);

	void Rotate(const Orientation &angle);
	void Move(float dx, float dy, float dz);
	void Fix();

	void Concate(const std::vector<Particle> &parts);

	/**
	 * @brief GetRotationRadius
	 * @return The distance from beginning of the center of coordinate system
	 * to the farthest point of particle.
	 */
	double GetRotationRadius() const;

	const complex &GetRefractiveIndex() const;
	void SetRefractiveIndex(const complex &value);

	const Angle3d &GetSymmetry() const;
	virtual void GetParticalFacetIdRange(Facet */*id*/, int &/*begin*/, int &/*end*/) const {}

	bool IsNonConvex() const;
	void Output();

public:
	bool isAggregated = false;
	Orientation rotAngle;

protected:
	Angle3d m_symmetry;		///< angle of particle symmetry

	complex m_refractiveIndex;	///< complex value of refractive index of the particle
	bool m_isNonConvex;

protected:
	void SetDefaultNormals();
	void SetDefaultCenters();
	void Reset();
	void SetSymmetry(double beta, double gamma, double alpha = 0);
	virtual void SetFacetParams() {}

private:
	void SetDParams();
	void RotateNormals();
	void SetFacetIndices();

private:
	LocalRotator m_rotator;
};

