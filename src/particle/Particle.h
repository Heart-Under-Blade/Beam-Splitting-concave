#pragma once

#include <math.h>
#include <vector>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"
#include "Rotator.h"

struct ParticleFacet
{
	Facet original; // facet with origin coordinates of points
	Facet actual; // facet with actual coordinates of points
};

/**
 * @brief The Particle class is the base class inherited
 * by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle : public Array<ParticleFacet>
{
public:
	Particle();
	Particle(int nFacets, const complex &refrIndex, bool isNonConvex = false);

	Facet *GetActualFacet(int i);
	void SetFromFile(const std::string &filename, double reduceSize = -1);

	void Rotate(const Orientation &orientation);
	void Move(float dx, float dy, float dz);
	void Scale(double ratio);

	void Concate(const std::vector<Particle> &parts);
	void RemoveFacet(int index);
	void CommitState();

	/**
	 * @brief LongRadius
	 * @return The distance from beginning of the center of coordinate system
	 * to the farthest point of particle.
	 */
	double LongRadius() const;

	/**
	 * @brief Area of the particle. Computes with summarising areas
	 * of all facets of the particle
	 * @return value of area
	 */
	double Area() const;

	/**
	 * @brief Compute the volume of the particle
	 * with splitting into tetrahedrons
	 * @return volume of the particle
	 */
	double Volume() const;

	/**
	 * @brief Distance between two the most distant vertices of the particle
	 * @return value of distance
	 */
	double MaximalDimension() const;

	/**
	 * @brief Geomertical center of the particle
	 * @return coordinates of center
	 */
	Point3f Center() const;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	/**
	 * @brief A sum of areas of each facet of the particle
	 * @return value of area
	 */
	double Area();

	double MaximalDimention() const;

	const complex &GetRefractiveIndex() const;
	void SetRefractiveIndex(const complex &value);

<<<<<<< HEAD
	const Symmetry &GetSymmetry() const;
	virtual void GetParticalFacetIdRangeByFacetId(
			int /*id*/, int &/*begin*/, int &/*end*/) const {}

	bool IsConcave() const;
=======
	const Orientation &GetSymmetry() const;
	virtual void GetParticalFacetIdRange(Facet */*id*/,
										 int &/*begin*/, int &/*end*/) const {}
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	bool IsNonConvex() const;
	void Output();

public:
	bool isAggregated = false;
	Orientation rotAngle;

protected:
	Orientation m_symmetry;		///< angle of particle symmetry

	// REF move this to somewhere
	complex m_refractiveIndex;	///< complex value of refractive index of the particle

	bool m_isNonConvex;

protected:
	void ReduceSmallEdges(double minSize);
	void SetDefaultNormals();
	void SetDefaultCenters();
	void ResetPosition();
	void SetSymmetry(double beta, double gamma);
	virtual void SetFacetParams() {}

private:
	void SetDParams();
	void RotateNormals();
	void SetFacetIndices();
	int ReduceEdge(int facetNo, int i1, int i2);

private:
	LocalRotator m_rotator;
};

