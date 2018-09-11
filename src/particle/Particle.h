#pragma once

#include <math.h>

#include "compl.hpp"
#include "geometry_lib.h"
#include "Facet.h"
#include <vector>

#define ROT_MTR_RANK 3

class Angle
{
public:
	double alpha;
	double beta;
	double gamma;

	Angle() {}
	Angle(double a, double b, double g)
	{
		alpha = a;
		beta = b;
		gamma = g;
	}

	void ToRadian()
	{
		alpha = DegToRad(alpha);
		beta = DegToRad(beta);
		gamma = DegToRad(gamma);
	}

	void ToDegree()
	{
		alpha = RadToDeg(alpha);
		beta = RadToDeg(beta);
		gamma = RadToDeg(gamma);
	}

	static double DegToRad(double deg)
	{
		return (deg*M_PI)/180;
	}

	static double RadToDeg(double rad)
	{
		return (rad*180)/M_PI;
	}

};

struct ParticleFacet
{
	Facet origin; // facet with origin coordinates of points
	Facet actual; // facet with actual coordinates of points
};

/**
 * @brief The Particle class
 * The base class inherited by other concrete particle classes.
 * Vertices are ordered by counterclock-wise direction if you see from outside.
 */
class Particle : public Array<ParticleFacet>
{
public:
	Particle();
	Particle(size_t nFacets, const complex &refrIndex, bool isNonConvex = false);

	Facet *GetActualFacet(size_t i);
	void SetFromFile(const std::string &filename);

	void Rotate(const Angle &angle);
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

	const Symmetry &GetSymmetry() const;
	virtual void GetParticalFacetIdRangeByFacetId(int /*id*/, int &/*begin*/, int &/*end*/) const {}

	bool IsNonConvex() const;
	void Output();

public:
	bool isAggregated = false;
	Angle rotAngle;

protected:
	Symmetry m_symmetry;		///< angle of particle symmetry

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
	void RotatePoint(const Point3f &point, Point3f &result);
	void SetRotateMatrix();
	void ReadSymmetry(const int bufSize, char *trash, char *buff,
					  std::ifstream pfile, char *ptr);
	void SetFacetIndices();

private:
	double m_rotMatrix[ROT_MTR_RANK][ROT_MTR_RANK];	///< rotation matrix for vertices
};

