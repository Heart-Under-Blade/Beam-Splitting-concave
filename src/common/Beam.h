#pragma once

#include "common.h"
#include "compl.hpp"
#include "JonesMatrix.h"
#include "float.h"
#include "BigInteger.hh"
#include "geometry_lib.h"
#include "Facet.h"
#include "Tracks.h"
#include "TrackTree.h"

template <class T>
class SplittedBeams
{
public:
	T internal;
	T external;

	void SetBeams(const Polygon &beamShape)
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
	bool IsInsideAtAct(int act) const
	{
		int mask = 1;
		mask <<= act;
		return !(locations & mask);
	}

	void RefreshId(int facetId, int nFacets)
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
//	TrackNode *node;

#ifdef MODE_FIXED_OR
	std::vector<Point3f> dirs;
	std::vector<Polygon> pols;
#endif
};

/**
 * @brief A plane-parallel optical beam that is created by
 * act of reflection / refraction when a light incidents on a Particle.
 */
struct Beam : public Polygon, public Track
{
public:
	Beam();
	Beam(const Beam &other);
	Beam(const Polygon &other);
	Beam(Beam &&other);

	Vector3f RotateSpherical(const Vector3f &dir, const Vector3f &polarBasis);

	void SetPolygon(const Polygon &other); // REF: мб просто искользовать "="?
	virtual void SetDefault();
	void CopyTrack(const Track &other);

	Beam & operator = (const Beam &other);
	Beam & operator = (const Polygon &other);
	Beam & operator = (Beam &&other);

	void SetIsInside(bool isIn);

	/**
	 * @brief Checks wheather the Beam is created by shadow of the Particle
	 * @return true if beam is shadow otherwise false
	 */
	bool IsShadow() const;

	void MultiplyByFresnel(const complex &f1, const complex &f2);
	void RotateJones(const Vector3f &normal, bool isInside);

	/**
	 * @brief Outputs beam params. Use it with std::cout or std::ofstream
	 * @param os output stream
	 * @param beam this beam
	 * @return output stream
	 */
	friend std::ostream & operator << (std::ostream &os, const Beam &beam);

public:
	Matrix2x2c Jones;	///< State of electric field of the Beam implemented as Jones's matrix
	Point3f direction;	///< Direction of propagation of the Beam
	Point3f polarizationBasis;	///< Basis of polarization vector

private:
	void SetDefault(Beam &other);
	void Copy(const Beam &other);
};

struct Axes
{
	Point3d horizontal;
	Point3d vertical;
};

struct VertexOrder
{
	int begin;
	int startIndex;
	int endIndex;
	int inc;

	VertexOrder() {}

	void SetOrder(bool isCcw, int nVertices)
	{
		if (isCcw)
		{
			begin = 0;
			startIndex = nVertices-1;
			endIndex = -1;
			inc = -1;
		}
		else
		{
			begin = nVertices-1;
			startIndex = 0;
			endIndex = nVertices;
			inc = 1;
		}
	}
};

struct BeamInfo
{
	bool isShadow;
	double area;
	double projLenght;
	double opticalLengths[3];
	Point3f beamBasis;
	Point3f centerf;
	Point3f normal;
	Point3d center;
	Point3d projectedCenter;
	Point3d normald;
	Point3d lenIndices;
	Beam *beam;
	Axes csAxes; // coordinate system axes
	VertexOrder order; // order of vertices
	std::vector<int> track;
};
