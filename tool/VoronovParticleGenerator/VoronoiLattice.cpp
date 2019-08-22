#include "VoronoiLattice.h"

#include <math.h>
#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>
#include <algorithm>

#include "geometry_lib.h"
#include "Converter.h"

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

VoronoiLattice::VoronoiLattice(double latticeSize, int splitRatio)
	: m_size(latticeSize)
{
	double r = latticeSize/2 - 40;
	int nRemovedLayers = 1;

	Lattice lattice(latticeSize, splitRatio, nRemovedLayers);

#ifdef _DEBUG // DEB
	std::ofstream dfile("debug.dat", std::ios::out);
#endif

	m_maxSiteDistance = 2.0*sqrt(2)*lattice.cellSize;
	int count = 0;

	int startLayer = nRemovedLayers;
	int endLayer = lattice.sideSize-nRemovedLayers;

	for (int i = startLayer; i < endLayer; ++i)
	{
		for (int j = startLayer; j < endLayer; ++j)
		{
			for (int k = startLayer; k < endLayer; ++k)
			{
				std::cout << "Progress: " << count++
						  << "/" << lattice.size << std::endl;

				Cell &cell = lattice.cells[i][j][k]; // lattice cell consists of facets

				std::vector<OrthoPlane> sitePlanes;
				DefineOrthogonalPlanes(cell, lattice, sitePlanes);

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

#ifdef _DEBUG // DEB
				double kk = 10;
				for (auto &p : sitePlanes)
				{
					dfile << p.center << std::endl;
		//			dfile << p.center + p.normal * kk << std::endl << std::endl;
				}
				dfile << std::endl;
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
	}

	file.close();
}

bool VoronoiLattice::RemoveDistantPlanes(const std::list<Point3f> &points,
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
		std::list<Point3f> &points)
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
#ifdef _DEBUG // DEB
					if (std::isnan(x.coordinates[0]) ||
							std::isnan(x.coordinates[1]) ||
							std::isnan(x.coordinates[2]))
						int d = 0;
#endif
					points.push_back(x);
				}
			}
		}
	}
}

void VoronoiLattice::RemoveExternalPoints(
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
	}

	if (points.size() < 3)
	{
		points.clear();
	}
}

void VoronoiLattice::RemoveDuplicatedPoints(std::list<Point3f> &points)
{
	for (auto it0 = points.begin(); it0 != points.end(); ++it0)
	{
		for (auto it1 = points.begin(); it1 != points.end(); ++it1)
		{
			if (it1 != it0)
			{
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

void VoronoiLattice::RemoveSameLinePoints(std::list<Point3f> &points)
{
	if (points.empty())
	{
		return;
	}

	std::vector<Point3f> _points{std::begin(points), std::end(points)};

	for (int i = 0; i < _points.size()-1; ++i)
	{
		for (int j = i + 1; j < _points.size(); ++j)
		{
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

	points = std::list<Point3f>{std::begin(_points), std::end(_points)};
}

void VoronoiLattice::OrderPoints(const Vector3f &planeNormal,
								 std::list<Point3f> &points) const
{
	if (!points.empty())
	{
		std::vector<Point3f> _points{std::begin(points), std::end(points)};
		Point3f pBase = _points.back();

		for (int iNext = 0; iNext < _points.size(); ++iNext)
		{
			Vector3f vBase = _points[iNext] - pBase;

			for (int i = iNext + 1; i <= _points.size(); ++i)
			{
				int iLast = (i == _points.size()) ? 0 : i;
				Vector3f vLast = _points[iLast] - pBase;
				Vector3f cp = CrossProduct(vBase, vLast);

				if (DotProduct(planeNormal, cp) < EPS_SAME_LINE)
				{
					std::swap(_points[iNext], _points[iLast]);
					vBase = _points[iNext] - pBase;
				}
			}

			pBase = _points[iNext];
		}

		points = std::list<Point3f>{std::begin(_points), std::end(_points)};
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

			bool isOk = true;
			Vector3f n = pBase.normal;
			n.coordinates[3] = -DotProduct(pBase.center, pBase.normal);
			Point3f x = ProjectPointToPlane(pCurr.center, v0ToEdge, n);

			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{pCurr.cell->id, edgeV, x});
#ifdef _DEBUG // DEB
				if (edgeLines.size() == 22)
					int d = 0;
#endif
			}
		}
	}
}

void VoronoiLattice::DefineOrthogonalPlanes(const Cell &baseCell,
											Lattice &lattice,
											std::vector<OrthoPlane> &sitePlanes)
{
	for (int i = 0; i < lattice.sideSize; ++i)
	{
		for (int j = 0; j < lattice.sideSize; ++j)
		{
			for (int k = 0; k < lattice.sideSize; ++k)
			{
				Cell &cell = lattice.cells[i][j][k];

				if (cell.id != baseCell.id)
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
					}
				}
			}
		}
	}

#ifdef _DEBUG // DEB
	int fff = 0;
#endif
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
