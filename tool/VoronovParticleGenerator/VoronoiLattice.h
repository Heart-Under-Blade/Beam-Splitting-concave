#pragma once

#include "geometry_lib.h"
#include <vector>
#include <list>

/**
 * @brief A plane between current site and other site
 */
struct OrthoPlane
{
	int id;
	Point3f normal;  ///< Normal of planes
	Point3f center;  ///< Center of planes
};

struct EdgeLine
{
	int planeId;
	Vector3f vector;
	Point3f beginPoint;
};

struct Intersection
{
	Point3f point;
	int firstPlaneId;
	int secondPlaneId;
};

typedef std::vector<std::list<Point3f>> Cell;

/**
 * @brief A 3D geometrical cubical lattice created by using of Voronoi
 * tesselation.
 * Consists of arbitrary formed particles (cells) connected each other with
 * common facets.
 */
class VoronoiLattice
{
public:
	/**
	 * @brief VoronoiLattice constructor. Generates a random Voronoi lattice
	 * with the central point in (0, 0, 0). Minimum point is
	 * (-latticeSize/2, -latticeSize/2, -latticeSize/2), maximum is
	 * (latticeSize/2, latticeSize/2, latticeSize/2).
	 * @param latticeSize size of the one dimention of the lattice in given
	 * units
	 * @param splitRatio number of divisions of the one dimention of lattice.
	 * The number of cells is splitRatio*splitRatio*splitRatio
	 */
	VoronoiLattice(double latticeSize, double splitRatio);

private:
	double m_size; ///< Lattice size

	/**
	 * @brief Generates sites of the lattice
	 * @param latticeSize size of the lattice in given units
	 * @param splitRatio
	 * @param sites generated sites
	 */
	void GenerateSites(double latticeSize, double splitRatio,
					   std::vector<Point3f> &sites);

	void DefineOrthogonalPlanes(int siteIndex, const std::vector<Point3f> &sites,
						   std::vector<OrthoPlane> &sitePlanes);

	void DefineFacetEdgeLines(const std::vector<OrthoPlane> &sitePlanes,
							  int planeIndex,
							  std::vector<EdgeLine> &edgeLines);

	void DefineIntersections(const std::vector<EdgeLine> &edgeLines,
							 const Vector3f &planeNormal,
							 std::list<Intersection> &points);

	/**
	 * @brief Remove intersection points are situated in space in front of each
	 * plane of the current site
	 * @param siteIndex index of the current site
	 * @param points intersection points that become vertices of the facet
	 */
	void RemoveExternalPoints(const std::vector<OrthoPlane> &sitePlanes,
							 std::list<Intersection> &points);

	void RemoveDuplicatedPoints(std::list<Intersection> &points);
	void RemoveSameLinePoints(std::list<Intersection> &points);

	void OrderPoints(const Vector3f &planeNormal,
					 std::list<Intersection> &points) const;

	void OutputFacets(std::vector<std::list<Point3f>> &facets);
	void FixPointsOrder(const Vector3f &normal, std::list<Point3f> &points);
	bool RemoveDistantPlanes(const std::list<Intersection> &points,
							const OrthoPlane &currentPlane,
							std::vector<OrthoPlane> &sitePlanes);
};
