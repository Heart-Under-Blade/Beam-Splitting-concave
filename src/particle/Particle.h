#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"
#include "Orientation.h"

#include <vector>

#define ROT_MTR_RANK 3

class Angle
{
public:
	double alpha;
	double beta;
	double gamma;
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
	Particle(int nFacets, bool isNonConvex = false);

	void SetFromFile(const std::string &filename);

	void Rotate(double beta, double gamma, double alpha);
	void Move(float dx, float dy, float dz);
	void Fix();
	void Resize(double size);
	void Concate(const std::vector<Particle> &parts);

	/**
	 * @brief LongRadius
	 * @return The distance from beginning of the center of coordinate system
	 * to the farthest point of particle.
	 */
	double LongRadius() const;

	double MaximalDimention() const;

	const complex &GetRefractiveIndex() const;
	void SetRefractiveIndex(const complex &value);

	virtual void GetParticalFacetIdRangeByFacetId(
			int /*id*/, int &/*begin*/, int &/*end*/) const {}

	bool IsConcave() const;

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

	const Orientation &GetSymmetry() const;
	virtual void GetPartByFacet(Facet */*facet*/, Array<Facet*> &facets);

	void Output();

public:
	Facet facets[MAX_FACET_NUM];	///< all facets of particle
	int nFacets;					///< number of facets
	bool isAggregated = false;

	Angle rotAngle;

protected:
	Facet defaultFacets[MAX_FACET_NUM];

	complex m_refractiveIndex;	///< complex value of refractive index of the particle
	Orientation m_symmetry;		///< angle of particle symmetry
	bool m_isNonConvex;

protected:
	void Init(int facetCount, const complex &refrIndex);

	void SetDefaultNormals();
	void SetDefaultCenters();

	void Reset();
	void Scale(double ratio);
	void ResetPosition();
	void SetSymmetry(double beta, double gamma);
	void GetFacets(int end, int begin, Array<Facet*> &facets);
	virtual void SetFacetParams() {}

private:
	void SetDParams();
	void RotateNormals();
	void RotatePoint(const Point3f &point, Point3f &result);
	void RotateCenters();
	void SetRotateMatrix(double beta, double gamma, double alpha);

private:
	double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
	void ReadSymmetry(const int bufSize, char *trash, char *buff,
					  std::ifstream pfile, char *ptr);
};

