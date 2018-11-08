#pragma once

#include "global.h"
#include "compl.hpp"
#include "JonesMatrix.h"
#include "float.h"
#include "BigInteger.hh"
#include "geometry_lib.h"
#include "Facet.h"
#include "Tracks.h"

template <class T>
class SplittedBeams
{
public:
	T internal;
	T external;

	void SetBeams(const Polygon1 &beamShape)
	{
		internal.Clear();
		internal.SetPolygon(beamShape);

		external.Clear();
		external.SetPolygon(beamShape);

#ifdef MODE_FIXED_OR
		internal.pols.push_back(beamShape);
		external.pols.push_back(beamShape);
#endif
	}
};

class Light
{
public:
	Point3f direction;
	Point3f polarizationBasis;
};

class Track
{
public:
	Track()
	{
		id = 0;
		actNo = -1;
		locations = 0;
	}

	void Update(Facet *face)
	{
		facet = face;
		++actNo;
	}

	/**
	 * @brief Shows where the Beam is located towards the Particle inside or outside
	 * @param iAct r/r act number
	 * @return true if beam location is "inside" on r/r act, otherwise returns false
	 */
	bool IsInsideOnAct(int iAct) const
	{
		int mask = 1;
		mask <<= iAct;
		return !(locations & mask);
	}

	void RecomputeTrackId(int facetId, int nFacets)
	{
		id = (id + (facetId + 1)) * (nFacets + 1);
	}

	Track & operator = (const Track &other)
	{
		if (this != &other) // OPT: попробовать убрать это уловие для ускорения
		{
			id = other.id;
			locations = other.locations;
			actNo = other.actNo;
			facet = other.facet;
#ifdef MODE_FIXED_OR
			pols = other.pols;
			dirs = other.dirs;
#endif
		}

		return *this;
	}

public:
	IdType id; ///< Unique id of beam calculated by facet ids
	int locations;	///< Every bit of the variable represents location of beam
					///< after an r/r act from the right to the left
					///< "0" when beam location is "inside" and "1" if it's "outside"

	int actNo; ///< Current r/r act number
	Facet *facet; ///< Last incident facet of the Particle

#ifdef MODE_FIXED_OR
	std::vector<Point3f> dirs;
	std::vector<Polygon1> pols;
#endif
};

/**
 * @brief A plane-parallel optical beam that is created by
 * act of reflection / refraction when a light incidents on a Particle.
 */
class Beam : public Polygon1, public Light, public Track
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(const Polygon1 &other);
	Beam(Beam &&other);

	Vector3f RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis);

	void SetPolygon(const Polygon1 &other); // REF: мб просто искользовать "="?
	virtual void Clear();
	void CopyTrack(const Track &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon1 &other);
	Beam & operator = (const Light &other);
	Beam & operator = (Beam &&other);

	void SetLocation(bool isIn);

	void MultiplyJonesMatrix(const complex &f1, const complex &f2);
	void RotateJones(const Vector3f &normal);

	// REF: перенести в Diffractor
	complex DiffractionIncline(const Point3d& pt, double wavelength) const; ///< calculate diffraction at the point /b pt
	//--------------------------

	void GetSpherical(double &fi, double &teta) const;

	/**
	 * @brief Outputs beam params. Use it with std::cout or std::ofstream
	 * @param os
	 * @param beam
	 * @return
	 */
	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

public:
	Matrix2x2c Jones;	///< Jones matrix of beam
	bool isInside; 		///< Beam state towards the particle (inside or outside)

private:
	void SetDefault(Beam &other);
	void Copy(const Beam &other);
};

