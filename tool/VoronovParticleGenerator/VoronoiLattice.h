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

typedef std::vector<std::vector<Point3f>> Matrix3x3p;
typedef std::vector<Matrix3x3p> Matrix3x3x3p;

struct SiteIndex
{
	int i;
	int j;
	int k;
};

class SiteLattice
{
public:
	int sideSize;
	int size;

	Matrix3x3x3p sites;

	const Point3f &GetSite(const SiteIndex &index) const
	{
		return sites[index.i][index.j][index.k];
	}

	SiteLattice(double latticeSize, int splitRatio)
	{
		sideSize = splitRatio;
		size = sideSize*sideSize*sideSize;

		double cellSize = latticeSize/sideSize;

		// add one layer of points to make borders
		latticeSize += cellSize * 2;
		sideSize += 2;

		double min = -latticeSize/2;
//		double max = latticeSize/2;

		Point3f from = Point3f(min, min, min);
		Point3f to = Point3f(min + cellSize,
							 min + cellSize,
							 min + cellSize);

		float &fromX = from.coordinates[0];
		float &fromY = from.coordinates[1];
		float &fromZ = from.coordinates[2];

		float &toX = to.coordinates[0];
		float &toY = to.coordinates[1];
		float &toZ = to.coordinates[2];

		for (int i = 0; i < sideSize; ++i)
		{
			sites.push_back(Matrix3x3p());
			double stepX = i*cellSize;

			for (int j = 0; j < sideSize; ++j)
			{
				sites[i].push_back(std::vector<Point3f>());
				double stepY = j*cellSize;

				for (int k = 0; k < sideSize; ++k)
				{
					double stepZ = k*cellSize;

					fromX += stepX;
					toX += stepX;

					fromY += stepY;
					toY += stepY;

					fromZ += stepZ;
					toZ += stepZ;

					Point3f p = RandomPoint(from, to);
					sites[i][j].push_back(p);
				}
			}
		}
	}

private:
	Point3f RandomPoint(const Point3f &from, const Point3f &to)
	{
		Point3f p;
		p.coordinates[0] = RandomValueF(from.coordinates[0], to.coordinates[0]);
		p.coordinates[1] = RandomValueF(from.coordinates[1], to.coordinates[1]);
		p.coordinates[2] = RandomValueF(from.coordinates[2], to.coordinates[2]);
		return p;
	}

	float RandomValueF(float from, float to)
	{
		return from + static_cast<float>(rand())/
				(static_cast<float>(RAND_MAX/(to-from)));
	}

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
	VoronoiLattice(double latticeSize, int splitRatio);

private:
	double m_size; ///< Lattice size

	/**
	 * @brief Generates sites of the lattice
	 * @param latticeSize size of the lattice in given units
	 * @param splitRatio
	 * @param sites generated sites
	 */
	void GenerateSites(double latticeSize, int splitRatio,
					   Matrix3x3x3p &sites);

	void DefineOrthogonalPlanes(const SiteIndex &siteIndex,
								const SiteLattice &sites,
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

	void OutputFacets(const std::string &filename,
					  std::vector<std::list<Point3f>> &facets);
	void FixPointsOrder(const Vector3f &normal, std::list<Point3f> &points);
	bool RemoveDistantPlanes(const std::list<Intersection> &points,
							const OrthoPlane &currentPlane,
							std::vector<OrthoPlane> &sitePlanes);
	void RemovePlane(int planeNo, std::vector<OrthoPlane> &sitePlanes);
};
