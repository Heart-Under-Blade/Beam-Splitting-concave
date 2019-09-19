#include "VoronoiLattice.h"

<<<<<<< HEAD
#include <math.h>
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>
#include <algorithm>

#include "geometry_lib.h"
#include "Converter.h"

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> origin/refactor
#define EPS_SAME_FACET 0.003
//#define EPS_SAME_FACET 0.43921
#define EPS_DUPLICATE 0.003
#define EPS_SAME_LINE 0.003
=======
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
// cell_2
//#define EPS_SAME_CELL 0.003
//#define EPS_SAME_FACET 0.0005

// cell_13
//#define EPS_SAME_CELL 0.1
//#define EPS_SAME_FACET 0.0005

#define EPS_SAME_CELL 0.0001
#define EPS_SAME_FACET 0.05
//#define EPS_SAME_FACET 0.0005
#define EPS_DUPLICATE 0.005
#define EPS_SAME_LINE 0.005

void VoronoiLattice::AddBorderPlane(std::vector<OrthoPlane> &sitePlanes,
									const Point3f &n, const Point3f &c,
									double r, Cell &cell)
{
	double len = Length(c - cell.site);

	if (len < m_maxSiteDistance)
	{
		sitePlanes.push_back(OrthoPlane{&cell, n, c, r});
	}
}
<<<<<<< HEAD
<<<<<<< HEAD
=======
#define EPS_SAME_FACET 0.003
//#define EPS_SAME_FACET 0.43921
#define EPS_DUPLICATE 0.003
#define EPS_SAME_LINE 0.003
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor

VoronoiLattice::VoronoiLattice(double latticeSize, int splitRatio)
	: m_size(latticeSize)
{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
	Lattice lattice(latticeSize, splitRatio);
=======
>>>>>>> origin/refactor
=======
	Lattice lattice(latticeSize, splitRatio);
=======
>>>>>>> origin/refactor
	double r = latticeSize/2 - 40;
	int nRemovedLayers = 1;

	Lattice lattice(latticeSize, splitRatio, nRemovedLayers);
<<<<<<< HEAD
<<<<<<< HEAD
=======
	Lattice lattice(latticeSize, splitRatio);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor

#ifdef _DEBUG // DEB
	std::ofstream dfile("debug.dat", std::ios::out);
#endif

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> origin/refactor
	int nRemovedLayers = 2;
	int count = 0;

	for (int i = nRemovedLayers; i < lattice.sideSize-nRemovedLayers; ++i)
	{
		for (int j = nRemovedLayers; j < lattice.sideSize-nRemovedLayers; ++j)
		{
			for (int k = nRemovedLayers; k < lattice.sideSize-nRemovedLayers; ++k)
=======
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
	m_maxSiteDistance = 2.0*sqrt(2)*lattice.cellSize;
	int count = 0;

	int startLayer = nRemovedLayers;
	int endLayer = lattice.sideSize-nRemovedLayers;

	for (int i = startLayer; i < endLayer; ++i)
	{
		for (int j = startLayer; j < endLayer; ++j)
		{
			for (int k = startLayer; k < endLayer; ++k)
<<<<<<< HEAD
<<<<<<< HEAD
=======
	int nRemovedLayers = 2;
	int count = 0;

	for (int i = nRemovedLayers; i < lattice.sideSize-nRemovedLayers; ++i)
	{
		for (int j = nRemovedLayers; j < lattice.sideSize-nRemovedLayers; ++j)
		{
			for (int k = nRemovedLayers; k < lattice.sideSize-nRemovedLayers; ++k)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
			{
				std::cout << "Progress: " << count++
						  << "/" << lattice.size << std::endl;

				Cell &cell = lattice.cells[i][j][k]; // lattice cell consists of facets

				std::vector<OrthoPlane> sitePlanes;
				DefineOrthogonalPlanes(cell, lattice, sitePlanes);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> origin/refactor
=======
=======
>>>>>>> origin/refactor
//				AddBorderPlane(sitePlanes, Point3f(0,0,1, r),
//							   Point3f(cell.site.coordinates[0],cell.site.coordinates[1],r), r, cell);
//				AddBorderPlane(sitePlanes, Point3f(0,1,0, r),
//							   Point3f(cell.site.coordinates[0],r,cell.site.coordinates[2]), r, cell);
//				AddBorderPlane(sitePlanes, Point3f(1,0,0, r),
//							   Point3f(r,cell.site.coordinates[1],cell.site.coordinates[2]), r, cell);
//				AddBorderPlane(sitePlanes, Point3f(0,0,-1, r),
//							   Point3f(cell.site.coordinates[0],cell.site.coordinates[1],-r), r, cell);
//				AddBorderPlane(sitePlanes, Point3f(0,-1,0, r),
//							   Point3f(cell.site.coordinates[0],-r,cell.site.coordinates[2]), r, cell);
//				AddBorderPlane(sitePlanes, Point3f(-1,0,0, r),
//							   Point3f(-r,cell.site.coordinates[1],cell.site.coordinates[2]), r, cell);

<<<<<<< HEAD
<<<<<<< HEAD
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
#ifdef _DEBUG // DEB
				double kk = 10;
				for (auto &p : sitePlanes)
				{
					dfile << p.center << std::endl;
		//			dfile << p.center + p.normal * kk << std::endl << std::endl;
				}
				dfile << std::endl;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> origin/refactor
#endif
				int planeCount = 0;

				for (int planeNo = 0; planeNo < sitePlanes.size(); ++planeNo)
				{
					std::vector<EdgeLine> edgeLines;
					DefineFacetEdgeLines(sitePlanes, planeNo, edgeLines);

#ifdef _DEBUG // DEB
					std::ofstream ffile("facet.dat", std::ios::out);
//					for (auto &pp : edgeLines)
//					{
//						Point3f b = (pp.vector * k) + pp.beginPoint;
//						ffile << b << std::endl;
//						ffile << (pp.vector * -k) + pp.beginPoint << std::endl << std::endl;
//					}
//					dfile << site.point << std::endl;
//					dfile << plane.center << std::endl;
					ffile << std::endl;
#endif
					auto &normal = sitePlanes[planeNo].normal;

					std::list<Intersection> points;
					DefineIntersections(edgeLines, normal, points);

#ifdef _DEBUG // DEB
					ffile << std::endl;
//					for (auto &pp : points)
//					{
//						ffile << pp.point << std::endl;
//					}
#endif
					RemoveExternalPoints(sitePlanes, points);
					RemoveDuplicatedPoints(points);
					RemoveSameLinePoints(points);
#ifdef _DEBUG // DEB
					if (planeNo == 3)
						int fff = 0;
#endif
					OrderPoints(normal, points);
		//			FixPointsOrder(normal, points);

		//			RemoveDistantPlanes(points, sitePlanes[planeId], sitePlanes);

					if (!points.empty())
					{
//						dfile << sitePlanes[planeNo].center << std::endl;

						std::list<Point3f> facet;

						for (auto &p : points)
						{
//							facet.push_back(p.point);
							facet.push_back(p);
						}

						cell.facets.push_back(facet);
						facet.reverse();
						sitePlanes[planeNo].cell->facets.push_back(facet);
					}
					else
					{
						RemovePlane(planeNo, sitePlanes);
						--planeNo;
					}

					++planeCount;
//					ffile.close();
//					std::cout << "count = " << std::to_string(count) << ", " << planeCount << std::endl;
=======
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
				if (count == 13)
				{
					for (auto &p : sitePlanes)
					{
						if (p.cell->id == 218)
							int fff = 0;
					}
				}
#endif
				int planeCount = 0;
				const auto &ids = cell.checkedSiteIds;

				for (int planeNo = 0; planeNo < sitePlanes.size(); ++planeNo)
				{
					auto &plane = sitePlanes[planeNo];
					auto it = std::find(std::begin(ids), std::end(ids),
										plane.cell->id);

					if (it == std::end(ids))
					{
						std::vector<EdgeLine> edgeLines;
						DefineFacetEdgeLines(sitePlanes, planeNo, edgeLines);
#ifdef _DEBUG // DEB
//						dfile << site.point << std::endl;
//						dfile << plane.center << std::endl;
						if (count == 13 && plane.cell->id == 218)
							int fff = 0;
#endif
						auto &normal = plane.normal;

						std::list<Point3f> points;
						DefineIntersections(edgeLines, normal, points);
#ifdef _DEBUG // DEB
						if (count == 13 && planeCount == 4)
						{
//							for (auto &pp : edgeLines)
//							{
//								Point3f b = (pp.vector * k) + pp.beginPoint;
//								ffile << b << std::endl;
//								ffile << (pp.vector * -k) + pp.beginPoint
//									  << std::endl << std::endl;
//							}
//							ffile << std::endl;
						}

						std::ofstream ffile("facet.dat", std::ios::out);
						if (plane.cell->id == 218)
						{
							for (auto &pp : points)
							{
								if (pp.coordinates[2] > 0 && pp.coordinates[2] < 120 &&
										pp.coordinates[1] > 0 && pp.coordinates[1] < 80 &&
										pp.coordinates[0] > 0 && pp.coordinates[0] < 120)
								{
									ffile << pp << std::endl;
								}
							}
						}
						ffile.close();
#endif
	//					std::cout << points.size() << ' ';
						RemoveExternalPoints(sitePlanes, points);
						RemoveDuplicatedPoints(points);

	//					std::cout << points.size() << ' ' << std::endl;
						RemoveSameLinePoints(points);
						OrderPoints(normal, points);
			//			FixPointsOrder(normal, points);
			//			RemoveDistantPlanes(points, sitePlanes[planeId], sitePlanes);

						if (!points.empty())
						{
	//						dfile << sitePlanes[planeNo].center << std::endl;
#ifdef _DEBUG // DEB
//							if (count == 13)
//								std::cout << plane.cell->id << std::endl;

//							if (cell.facets.size() == 8)
//								int ggg = 0;
#endif
							std::list<Point3f> facet;

							for (auto &p : points)
							{
								facet.push_back(p);
							}

							cell.facets.push_back(facet);
							facet.reverse();
							plane.cell->facets.push_back(facet);
						}
						else
						{
//							RemovePlane(planeNo, sitePlanes);
//							--planeNo;
						}

						++planeCount;
	//					std::cout << "count = " << std::to_string(count) << ", " << planeCount << std::endl;
					}
<<<<<<< HEAD
<<<<<<< HEAD
				}

				auto facets = lattice.ToFacets(&cell);
				std::string filename = "cell_" + std::to_string(count);
				OutputFacets(filename + ".dat", facets, false);
#ifdef _DEBUG // DEB
				if (count == 19)
					int uuu = 0;
#endif
				Converter converter;
				std::vector<Facet> triangles;
				converter.Triangulate(facets, triangles);
=======
#endif
				int planeCount = 0;

				for (int planeNo = 0; planeNo < sitePlanes.size(); ++planeNo)
				{
					std::vector<EdgeLine> edgeLines;
					DefineFacetEdgeLines(sitePlanes, planeNo, edgeLines);

#ifdef _DEBUG // DEB
					std::ofstream ffile("facet.dat", std::ios::out);
//					for (auto &pp : edgeLines)
//					{
//						Point3f b = (pp.vector * k) + pp.beginPoint;
//						ffile << b << std::endl;
//						ffile << (pp.vector * -k) + pp.beginPoint << std::endl << std::endl;
//					}
//					dfile << site.point << std::endl;
//					dfile << plane.center << std::endl;
					ffile << std::endl;
#endif
					auto &normal = sitePlanes[planeNo].normal;

					std::list<Intersection> points;
					DefineIntersections(edgeLines, normal, points);

#ifdef _DEBUG // DEB
					ffile << std::endl;
//					for (auto &pp : points)
//					{
//						ffile << pp.point << std::endl;
//					}
#endif
					RemoveExternalPoints(sitePlanes, points);
					RemoveDuplicatedPoints(points);
					RemoveSameLinePoints(points);
#ifdef _DEBUG // DEB
					if (planeNo == 3)
						int fff = 0;
#endif
					OrderPoints(normal, points);
		//			FixPointsOrder(normal, points);

		//			RemoveDistantPlanes(points, sitePlanes[planeId], sitePlanes);

					if (!points.empty())
					{
//						dfile << sitePlanes[planeNo].center << std::endl;

						std::list<Point3f> facet;

						for (auto &p : points)
						{
//							facet.push_back(p.point);
							facet.push_back(p);
						}

						cell.facets.push_back(facet);
						facet.reverse();
						sitePlanes[planeNo].cell->facets.push_back(facet);
					}
					else
					{
						RemovePlane(planeNo, sitePlanes);
						--planeNo;
					}

					++planeCount;
//					ffile.close();
//					std::cout << "count = " << std::to_string(count) << ", " << planeCount << std::endl;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
				}

#ifdef _DEBUG // DEB
//				dfile.close();
#endif
				auto facets = lattice.ToFacets(&cell);
				std::string filename = "cell_" + std::to_string(count);
				OutputFacets(filename + ".dat", facets);
				Converter converter;
				std::vector<Facet> triangles;
				converter.Triangulate(lattice.ToFacets(&cell), triangles);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
				converter.WriteStl(triangles, filename);
			}
		}
	}

//	for (auto &cell : lattice)
//	{
//		OutputFacets("cell_" + std::to_string(count), cell);
//	}

	Converter converter;
	std::vector<Facet> triangles;

	std::vector<Facet> facet;
	converter.Triangulate(lattice.ToFacets(), facet);
	triangles.insert(triangles.end(), facet.begin(), facet.end());

	converter.WriteStl(triangles, "voronoi");
}

void VoronoiLattice::RemovePlane(int planeNo, std::vector<OrthoPlane> &sitePlanes)
{
	auto it = sitePlanes.begin();
	std::advance(it, planeNo);
	sitePlanes.erase(it);
}

void VoronoiLattice::OutputFacets(const std::string &filename,
<<<<<<< HEAD
								  std::vector<std::vector<Point3f>> &facets,
								  bool isForGrapher)
{
	std::ofstream file(filename, std::ios::out);

	if (isForGrapher)
	{
		// size pf the maximal vector
		int maxSize = 0;
		for (int i = 0; i < facets.size(); ++i)
		{
			int newSize = facets[i].size();

			if (newSize > maxSize)
			{
				maxSize = newSize;
			}
		}

		for (int i = 0; i <= maxSize; ++i)
		{
			for (int j = 0; j < facets.size(); ++j)
			{
				if (i > facets[j].size())
				{
					file << "# # # # ";
				}
				else if (i == facets[j].size())
				{
					file << facets[j][0] << " # ";
				}
				else
				{
					Point3f &p = facets[j][i];
					file << p << " # ";
				}
			}

			file << std::endl;
		}
	}
	else
	{
		Point3f center = Point3f(0, 0, 0);
		int count = 0;

		for (int i = 0; i < facets.size(); ++i)
		{
			for (int j = 0; j < facets[i].size(); ++j)
			{
				center = center + facets[i][j];
				++count;
			}
		}

		center = center/count;

		for (int i = 0; i < facets.size(); ++i)
		{
			FixPointsOrder(facets[i][0] - center, facets[i]);

			for (int j = 0; j < facets[i].size(); ++j)
			{
				file << facets[i][j] << " " << std::endl;
			}

			file << std::endl;
		}
=======
								  std::vector<std::vector<Point3f>> &facets)
{
	// size pf the maximal vector
	int maxSize = 0;
	for (int i = 0; i < facets.size(); ++i)
	{
		int newSize = facets[i].size();

		if (newSize > maxSize)
		{
			maxSize = newSize;
		}
	}

	std::ofstream file(filename, std::ios::out);

	for (int i = 0; i <= maxSize; ++i)
	{
		for (int j = 0; j < facets.size(); ++j)
		{
			if (i > facets[j].size())
			{
				file << "# # # # ";
			}
			else if (i == facets[j].size())
			{
				file << facets[j][0] << " # ";
			}
			else
			{
				Point3f &p = facets[j][i];
				file << p << " # ";
			}
		}

		file << std::endl;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	}

	file.close();
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
bool VoronoiLattice::RemoveDistantPlanes(const std::list<Point3f> &points,
=======
bool VoronoiLattice::RemoveDistantPlanes(const std::list<Intersection> &points,
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
bool VoronoiLattice::RemoveDistantPlanes(const std::list<Intersection> &points,
=======
bool VoronoiLattice::RemoveDistantPlanes(const std::list<Point3f> &points,
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
										 const OrthoPlane &currentPlane,
										 std::vector<OrthoPlane> &sitePlanes)
{
	int originalSize = sitePlanes.size();

	if (originalSize > 0)
	{
		int count = 0;

		auto isDistant = [&](OrthoPlane &plane)
		{
			auto vToPlane = plane.center - currentPlane.center;
			bool flag = DotProduct(currentPlane.normal, vToPlane) > 0
					/*&& DotProduct(currentPlane.normal, plane.normal) > 0*/;

			if (flag)
			{
				for (const auto &p : points)
				{
//					if (count == p.firstPlaneId || count == p.secondPlaneId)
					{	// plane is intersect ortho plane
						flag = false;
						break;
					}
				}
			}

			++count;
			return flag;
		};

		sitePlanes.erase(std::remove_if(std::begin(sitePlanes),
										std::end(sitePlanes),
										isDistant),
						 sitePlanes.end());
	}

	bool isRemoved = sitePlanes.size() < originalSize;
	return isRemoved;
}

void VoronoiLattice::DefineIntersections(
		const std::vector<EdgeLine> &edgeLines, const Vector3f &planeNormal,
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		std::list<Point3f> &points)
=======
		std::list<Intersection> &points)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
		std::list<Intersection> &points)
=======
		std::list<Point3f> &points)
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
{
	for (int i = 0; i < edgeLines.size(); ++i)
	{
		for (int j = i + 1; j < edgeLines.size(); ++j)
		{
			if (j != i)
			{
				bool isOk = false;
				Point3f x = IntersectVectors(
							edgeLines[i].beginPoint, edgeLines[j].beginPoint,
							edgeLines[i].vector, edgeLines[j].vector,
							planeNormal, isOk);
				if (isOk)
				{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//					points.push_back(Intersection{x, 0, 0});
=======
>>>>>>> origin/refactor
=======
//					points.push_back(Intersection{x, 0, 0});
=======
>>>>>>> origin/refactor
#ifdef _DEBUG // DEB
					if (std::isnan(x.coordinates[0]) ||
							std::isnan(x.coordinates[1]) ||
							std::isnan(x.coordinates[2]))
						int d = 0;
#endif
<<<<<<< HEAD
<<<<<<< HEAD
=======
//					points.push_back(Intersection{x, 0, 0});
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
					points.push_back(x);
				}
			}
		}
	}
}

void VoronoiLattice::RemoveExternalPoints(
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
=======
>>>>>>> origin/refactor
		const std::vector<OrthoPlane> &sitePlanes,
		std::list<Intersection> &points)
{
	for (int i = 0; i < sitePlanes.size() && !points.empty(); ++i)
	{
		const Point3f &n = sitePlanes[i].normal;
		const Point3f &c = sitePlanes[i].center;
		points.remove_if([&](Intersection &p) {
//			Vector3f v = p.point-c;
			Vector3f v = p-c;
//			Normalize(v);
			double vv = DotProduct(n, v);
			bool res = vv > EPS_SAME_FACET;
			return res;
//			return DotProduct(n, p.point-c) > EPS_SAME_FACET;
		});
=======
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		const std::vector<OrthoPlane> &sitePlanes, std::list<Point3f> &points)
{
	for (auto &plane : sitePlanes)
	{
		Point3f n = plane.normal;
		n.coordinates[3] = plane.dParam;

//		double dd = plane.normal.coordinates[3];
		if (points.size() >= 3)
		{
			points.remove_if([&](Point3f &p) {
				Vector3f v = p - plane.center;
				double cosA = DotProduct(n, v);

				bool res = false;

				if (cosA > EPS_SAME_CELL)
				{
					Point3f x = ProjectPointToPlane(p, -n, n);
					double h = Length(x - p);
					res = h > EPS_SAME_FACET;
				}

				return res;
			});
		}
		else
		{
			break;
		}
<<<<<<< HEAD
<<<<<<< HEAD
=======
		const std::vector<OrthoPlane> &sitePlanes,
		std::list<Intersection> &points)
{
	for (int i = 0; i < sitePlanes.size() && !points.empty(); ++i)
	{
		const Point3f &n = sitePlanes[i].normal;
		const Point3f &c = sitePlanes[i].center;
		points.remove_if([&](Intersection &p) {
//			Vector3f v = p.point-c;
			Vector3f v = p-c;
//			Normalize(v);
			double vv = DotProduct(n, v);
			bool res = vv > EPS_SAME_FACET;
			return res;
//			return DotProduct(n, p.point-c) > EPS_SAME_FACET;
		});
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
	}

	if (points.size() < 3)
	{
		points.clear();
	}
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
void VoronoiLattice::RemoveDuplicatedPoints(std::list<Point3f> &points)
=======
void VoronoiLattice::RemoveDuplicatedPoints(std::list<Intersection> &points)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
void VoronoiLattice::RemoveDuplicatedPoints(std::list<Intersection> &points)
=======
void VoronoiLattice::RemoveDuplicatedPoints(std::list<Point3f> &points)
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
{
	for (auto it0 = points.begin(); it0 != points.end(); ++it0)
	{
		for (auto it1 = points.begin(); it1 != points.end(); ++it1)
		{
			if (it1 != it0)
			{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//				double len = Length((*it1).point - (*it0).point);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//				double len = Length((*it1).point - (*it0).point);
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//				double len = Length((*it1).point - (*it0).point);
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
				double len = Length((*it1) - (*it0));

				if (len < EPS_DUPLICATE)
				{
					points.erase(it1);
					it0 = points.begin();
					break;
				}
			}
		}
	}

	if (points.size() < 3)
	{
		points.clear();
	}
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
void VoronoiLattice::RemoveSameLinePoints(std::list<Point3f> &points)
=======
void VoronoiLattice::RemoveSameLinePoints(std::list<Intersection> &points)
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
void VoronoiLattice::RemoveSameLinePoints(std::list<Intersection> &points)
=======
void VoronoiLattice::RemoveSameLinePoints(std::list<Point3f> &points)
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
void VoronoiLattice::RemoveSameLinePoints(std::list<Intersection> &points)
=======
void VoronoiLattice::RemoveSameLinePoints(std::list<Point3f> &points)
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
{
	if (points.empty())
	{
		return;
	}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
	std::vector<Point3f> _points{std::begin(points), std::end(points)};
=======
	std::vector<Intersection> _points{std::begin(points), std::end(points)};
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
	std::vector<Intersection> _points{std::begin(points), std::end(points)};
=======
	std::vector<Point3f> _points{std::begin(points), std::end(points)};
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor

	for (int i = 0; i < _points.size()-1; ++i)
	{
		for (int j = i + 1; j < _points.size(); ++j)
		{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//			Vector3f vBase = _points[j].point - _points[i].point;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//			Vector3f vBase = _points[j].point - _points[i].point;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//			Vector3f vBase = _points[j].point - _points[i].point;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
			Vector3f vBase = _points[j] - _points[i];
			int k = j + 1;

			while (k <= _points.size())
			{
				if (k == _points.size())
				{
					if (i != 0)
					{
						k = 0;
					}
					else
					{
						break;
					}
				}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//				Vector3f vLast = _points[k].point - _points[i].point;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//				Vector3f vLast = _points[k].point - _points[i].point;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//				Vector3f vLast = _points[k].point - _points[i].point;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
				Vector3f vLast = _points[k] - _points[i];
				Vector3f cp = CrossProduct(vBase, vLast);
				bool isErased = Length(cp) < EPS_SAME_LINE;

				if (isErased)
				{
					if (Length(vBase) < Length(vLast))
					{
						auto itNext = _points.begin();
						std::advance(itNext, j);
						_points.erase(itNext);
						vBase = vLast;
					}
					else
					{
						auto itLast = _points.begin();
						std::advance(itLast, k);
						_points.erase(itLast);
					}
				}

				if (k == 0)
				{
					break;
				}

				if (!isErased)
				{
					++k;
				}
			}
		}
	}

	if (_points.size() < 3)
	{
		_points.clear();
	}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
	points = std::list<Point3f>{std::begin(_points), std::end(_points)};
=======
=======
>>>>>>> origin/refactor
	points = std::list<Intersection>{std::begin(_points), std::end(_points)};
>>>>>>> origin/refactor
}

void VoronoiLattice::OrderPoints(const Vector3f &planeNormal,
								 std::list<Point3f> &points) const
{
	if (!points.empty())
	{
<<<<<<< HEAD
		std::vector<Point3f> _points{std::begin(points), std::end(points)};
=======
	points = std::list<Intersection>{std::begin(_points), std::end(_points)};
=======
		std::vector<Intersection> _points{std::begin(points), std::end(points)};
//		Point3f pBase = _points.back().point;
=======
	points = std::list<Point3f>{std::begin(_points), std::end(_points)};
>>>>>>> origin/refactor
}

void VoronoiLattice::OrderPoints(const Vector3f &planeNormal,
								 std::list<Intersection> &points) const
{
	if (!points.empty())
	{
<<<<<<< HEAD
		std::vector<Intersection> _points{std::begin(points), std::end(points)};
//		Point3f pBase = _points.back().point;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
		std::vector<Point3f> _points{std::begin(points), std::end(points)};
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		Point3f pBase = _points.back();

		for (int iNext = 0; iNext < _points.size(); ++iNext)
		{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//			Vector3f vBase = _points[iNext].point - pBase;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//			Vector3f vBase = _points[iNext].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//			Vector3f vBase = _points[iNext].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
			Vector3f vBase = _points[iNext] - pBase;

			for (int i = iNext + 1; i <= _points.size(); ++i)
			{
				int iLast = (i == _points.size()) ? 0 : i;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//				Vector3f vLast = _points[iLast].point - pBase;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//				Vector3f vLast = _points[iLast].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//				Vector3f vLast = _points[iLast].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
				Vector3f vLast = _points[iLast] - pBase;
				Vector3f cp = CrossProduct(vBase, vLast);

				if (DotProduct(planeNormal, cp) < EPS_SAME_LINE)
				{
					std::swap(_points[iNext], _points[iLast]);
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
//					vBase = _points[iNext].point - pBase;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//					vBase = _points[iNext].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
=======
//					vBase = _points[iNext].point - pBase;
=======
>>>>>>> origin/feature/voronoi
>>>>>>> origin/refactor
					vBase = _points[iNext] - pBase;
				}
			}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			pBase = _points[iNext];
		}

		points = std::list<Point3f>{std::begin(_points), std::end(_points)};
=======
//			pBase = _points[iNext].point;
			pBase = _points[iNext];
		}

		points = std::list<Intersection>{std::begin(_points), std::end(_points)};
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
//			pBase = _points[iNext].point;
			pBase = _points[iNext];
		}

		points = std::list<Intersection>{std::begin(_points), std::end(_points)};
=======
			pBase = _points[iNext];
		}

		points = std::list<Point3f>{std::begin(_points), std::end(_points)};
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
	}
}

void VoronoiLattice::DefineFacetEdgeLines(
		const std::vector<OrthoPlane> &sitePlanes, int planeNo,
		std::vector<EdgeLine> &edgeLines)
{
	auto &pBase = sitePlanes[planeNo];

	for (int i = 0; i < sitePlanes.size(); ++i)
	{
		if (i != planeNo)
		{
			auto &pCurr = sitePlanes[i];
<<<<<<< HEAD
			Point3f edgeV = CrossProduct(pBase.normal, pCurr.normal);

			if (fabs(edgeV.coordinates[0]) +
					fabs(edgeV.coordinates[1]) +
					fabs(edgeV.coordinates[2]) < 3*FLT_EPSILON)
			{
				if (DotProduct(pBase.normal, pCurr.normal) < 0)
				{
					continue;
				}
				else
				{
					if (pBase.distance > pCurr.distance)
					{
						edgeLines.clear();
						return;
					}
					else
					{
						continue;
					}
				}
			}

			Point3f v0ToEdge = CrossProduct(pCurr.normal, edgeV);

=======

			Point3f edgeV = CrossProduct(pBase.normal, pCurr.normal);
<<<<<<< HEAD
<<<<<<< HEAD
=======
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor

			Point3f v0ToEdge = CrossProduct(pCurr.normal, edgeV);
//			Point3f v1ToEdge = CrossProduct(edgeV, pBase.normal);

//			bool isOk = false;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
=======
			Point3f v0ToEdge = CrossProduct(pCurr.normal, edgeV);

>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			bool isOk = true;
			Vector3f n = pBase.normal;
			n.coordinates[3] = -DotProduct(pBase.center, pBase.normal);
			Point3f x = ProjectPointToPlane(pCurr.center, v0ToEdge, n);
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD

			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{pCurr.cell->id, edgeV, x});
#ifdef _DEBUG // DEB
				if (edgeLines.size() == 22)
					int d = 0;
#endif
=======
=======
>>>>>>> origin/refactor
//			Point3f x = IntersectVectors(pCurr.center, pBase.center,
//										 v0ToEdge, v1ToEdge, edgeV, isOk);
			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{edgeV, x});
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
//			Point3f x = IntersectVectors(pCurr.center, pBase.center,
//										 v0ToEdge, v1ToEdge, edgeV, isOk);
			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{edgeV, x});
=======
>>>>>>> origin/refactor
=======

			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{pCurr.cell->id, edgeV, x});
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			}
		}
	}
}

void VoronoiLattice::DefineOrthogonalPlanes(const Cell &baseCell,
											Lattice &lattice,
											std::vector<OrthoPlane> &sitePlanes)
{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
	int maxSiteDistance = 2.0*sqrt(2)*lattice.cellSize;
	const auto &ids = baseCell.checkedSiteIds;

>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
	int maxSiteDistance = 2.0*sqrt(2)*lattice.cellSize;
	const auto &ids = baseCell.checkedSiteIds;

=======
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
	for (int i = 0; i < lattice.sideSize; ++i)
	{
		for (int j = 0; j < lattice.sideSize; ++j)
		{
			for (int k = 0; k < lattice.sideSize; ++k)
			{
				Cell &cell = lattice.cells[i][j][k];

				if (cell.id != baseCell.id)
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
				{	// site has not been checked yet
					Point3f normal = cell.site - baseCell.site;
					double len = Length(normal);

					if (len < m_maxSiteDistance)
					{
						Normalize(normal);

						// D-param
						Point3f center = (baseCell.site + cell.site)/2;
						normal.coordinates[3] = -DotProduct(center, normal);

						sitePlanes.push_back(OrthoPlane{&cell, normal, center,
														normal.coordinates[3],
											 Length(baseCell.site - cell.site)});
						cell.checkedSiteIds.push_back(baseCell.id);
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
				{
					auto it = std::find(std::begin(ids), std::end(ids), cell.id);

					if (it == std::end(ids))
					{	// site has not been checked yet
						Point3f normal = cell.site - baseCell.site;
						double len = Length(normal);

						if (len < maxSiteDistance)
						{
							Normalize(normal);

							// D-param
							Point3f center = (baseCell.site + cell.site)/2;
							normal.coordinates[3] = -DotProduct(center, normal);

							sitePlanes.push_back(OrthoPlane{&cell, normal, center});
							cell.checkedSiteIds.push_back(baseCell.id);
						}
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
=======
				{	// site has not been checked yet
					Point3f normal = cell.site - baseCell.site;
					double len = Length(normal);

					if (len < m_maxSiteDistance)
					{
						Normalize(normal);

						// D-param
						Point3f center = (baseCell.site + cell.site)/2;
						normal.coordinates[3] = -DotProduct(center, normal);

						sitePlanes.push_back(OrthoPlane{&cell, normal, center,
														normal.coordinates[3]});
						cell.checkedSiteIds.push_back(baseCell.id);
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
					}
				}
			}
		}
	}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
#ifdef _DEBUG // DEB
	int fff = 0;
#endif
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
=======
#ifdef _DEBUG // DEB
	int fff = 0;
#endif
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
	// sort by distance
	std::sort(sitePlanes.begin(), sitePlanes.end(),
			  [&](OrthoPlane &p1, OrthoPlane &p2) {
		return p1.normal.coordinates[3] < p2.normal.coordinates[3];
	});
}

void VoronoiLattice::FixPointsOrder(const Vector3f &normal,
									std::list<Point3f> &points)
{
	if (points.size() > 0)
	{
		Facet f;

		for (auto &p : points)
		{
			f.AddVertex(p);
		}

		Vector3f n = f.Normal();

		if (DotProduct(n, normal) < 0)
		{
			points.reverse();
		}
	}
}
<<<<<<< HEAD

void VoronoiLattice::FixPointsOrder(const Vector3f &normal,
									std::vector<Point3f> &points)
{
	if (points.size() > 0)
	{
		Facet f;

		for (auto &p : points)
		{
			f.AddVertex(p);
		}

		Vector3f n = f.Normal();

		if (DotProduct(n, normal) < 0)
		{
			for (int i = 0; i < points.size()/2; ++i)
			{
				Point3f buff = points[i];
				points[i] = points[points.size()-i];
				points[points.size()-i] = buff;
			}
		}
	}
}
=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
