#pragma once

#include "geometry_lib.h"
#include <vector>
#include <list>

struct EdgeLine
{
	int facetNo;
	Vector3f vector;
	Point3f beginPoint;
};

struct Cell
{
	int id;	///< id of site or cell
	Point3f site;	///< base point for generated point
	std::vector<std::list<Point3f>> facets; ///< facets of the Cell
	std::vector<int> checkedSiteIds;	///< ids of sites that facet had already
										/// generated for or imposible
										/// to generate one
};

typedef std::vector<std::vector<Cell>> Matrix3x3s;
typedef std::vector<Matrix3x3s> Matrix3x3x3s;

/**
 * @brief A plane between current site and other site
 */
struct OrthoPlane
{
	Cell *cell;
	Point3f normal;  ///< Normal of planes
	Point3f center;  ///< Center of planes
	double dParam;
};

struct SiteIndex
{
	int i;
	int j;
	int k;
};

class Lattice
{
public:
	int sideSize;
	int size;
	double cellSize;

	Matrix3x3x3s cells;

	std::vector<std::vector<Point3f>> ToFacets() const
	{
		std::vector<std::vector<Point3f>> facets;

		for (int i = 0; i < sideSize; ++i)
		{
			for (int j = 0; j < sideSize; ++j)
			{
				for (int k = 0; k < sideSize; ++k)
				{
					for (auto &f : cells[i][j][k].facets)
					{
						std::vector<Point3f> facet;

						for (auto &p : f)
						{
							facet.push_back(p);
						}

						facets.push_back(facet);
					}
				}
			}
		}

		return facets;
	}

	std::vector<std::vector<Point3f>> ToFacets(Cell* cell) const
	{
		std::vector<std::vector<Point3f>> facets;

		for (auto &f : cell->facets)
		{
			std::vector<Point3f> facet;

			for (auto &p : f)
			{
				facet.push_back(p);
			}

			facets.push_back(facet);
		}

		return facets;
	}

	const Cell &GetCell(const SiteIndex &index) const
	{
		return cells[index.i][index.j][index.k];
	}

	Lattice(double latticeSize, int splitRatio, int removedLayers)
	{
		int sr = splitRatio - removedLayers;
		size = sr*sr*sr;
		cellSize = latticeSize/splitRatio;

		// add one layer of points to make borders
		splitRatio += 2;
		double extendedLatticeSize = latticeSize + cellSize*2;

		sideSize = splitRatio;

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

		int count = 0;

		for (int i = 0; i < sideSize; ++i)
		{
			cells.push_back(Matrix3x3s());

			for (int j = 0; j < sideSize; ++j)
			{
				cells[i].push_back(std::vector<Cell>());

				for (int k = 0; k < sideSize; ++k)
				{
					fromZ += cellSize;
					toZ += cellSize;

					Cell cell;
					cell.id = count++;
					cell.site = RandomPoint(from, to);
					cells[i][j].push_back(cell);
#ifdef _DEBUG // DEB
					if (fabs(cell.site.coordinates[0]) > extendedLatticeSize ||
							fabs(cell.site.coordinates[1]) > extendedLatticeSize ||
							fabs(cell.site.coordinates[2]) > extendedLatticeSize)
						int fff = 0;
#endif
				}

				fromZ = min;
				toZ = min + cellSize;

				fromY += cellSize;
				toY += cellSize;
			}

			fromY = min;
			toY = min + cellSize;

			fromX += cellSize;
			toX += cellSize;
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
					   Matrix3x3x3s &sites);

	void DefineOrthogonalPlanes(const Cell &baseCell, Lattice &cells,
								std::vector<OrthoPlane> &sitePlanes);

	void DefineFacetEdgeLines(const std::vector<OrthoPlane> &sitePlanes,
							  int planeNo, std::vector<EdgeLine> &edgeLines);

	void DefineIntersections(const std::vector<EdgeLine> &edgeLines,
							 const Vector3f &planeNormal,
							 std::list<Point3f> &points);

	/**
	 * @brief Remove intersection points are situated in space in front of each
	 * plane of the current site
	 * @param siteIndex index of the current site
	 * @param points intersection points that become vertices of the facet
	 */
	void RemoveExternalPoints(const std::vector<OrthoPlane> &sitePlanes,
							 std::list<Point3f> &points);

	void RemoveDuplicatedPoints(std::list<Point3f> &points);
	void RemoveSameLinePoints(std::list<Point3f> &points);

	void OrderPoints(const Vector3f &planeNormal,
					 std::list<Point3f> &points) const;

	void OutputFacets(const std::string &filename,
					  std::vector<std::vector<Point3f>> &facets);
	void FixPointsOrder(const Vector3f &normal, std::list<Point3f> &points);
	bool RemoveDistantPlanes(const std::list<Point3f> &points,
							const OrthoPlane &currentPlane,
							std::vector<OrthoPlane> &sitePlanes);
	void RemovePlane(int planeNo, std::vector<OrthoPlane> &sitePlanes);
};
