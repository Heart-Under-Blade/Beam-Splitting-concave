#include "VoronoiLattice.h"

#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>
#include <algorithm>

#include "geometry_lib.h"
#include "Converter.h"

#define EPS_SAME_FACET 0.003
//#define EPS_SAME_FACET 0.43921
#define EPS_DUPLICATE 0.003
#define EPS_SAME_LINE 0.003

VoronoiLattice::VoronoiLattice(double latticeSize, int splitRatio)
	: m_size(latticeSize)
{
	Lattice lattice(latticeSize, splitRatio);

#ifdef _DEBUG // DEB
	std::ofstream dfile("debug.dat", std::ios::out);
#endif

	int nRemovedLayers = 2;
	int count = 0;

	for (int i = nRemovedLayers; i < lattice.sideSize-nRemovedLayers; ++i)
	{
		for (int j = nRemovedLayers; j < lattice.sideSize-nRemovedLayers; ++j)
		{
			for (int k = nRemovedLayers; k < lattice.sideSize-nRemovedLayers; ++k)
			{
				std::cout << "Progress: " << count++
						  << "/" << lattice.size << std::endl;

				Cell &cell = lattice.cells[i][j][k]; // lattice cell consists of facets

				std::vector<OrthoPlane> sitePlanes;
				DefineOrthogonalPlanes(cell, lattice, sitePlanes);

#ifdef _DEBUG // DEB
				double kk = 10;
				for (auto &p : sitePlanes)
				{
					dfile << p.center << std::endl;
		//			dfile << p.center + p.normal * kk << std::endl << std::endl;
				}
				dfile << std::endl;
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
	}

	file.close();
}

bool VoronoiLattice::RemoveDistantPlanes(const std::list<Intersection> &points,
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
		std::list<Intersection> &points)
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
//					points.push_back(Intersection{x, 0, 0});
					points.push_back(x);
				}
			}
		}
	}
}

void VoronoiLattice::RemoveExternalPoints(
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
	}

	if (points.size() < 3)
	{
		points.clear();
	}
}

void VoronoiLattice::RemoveDuplicatedPoints(std::list<Intersection> &points)
{
	for (auto it0 = points.begin(); it0 != points.end(); ++it0)
	{
		for (auto it1 = points.begin(); it1 != points.end(); ++it1)
		{
			if (it1 != it0)
			{
//				double len = Length((*it1).point - (*it0).point);
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

void VoronoiLattice::RemoveSameLinePoints(std::list<Intersection> &points)
{
	if (points.empty())
	{
		return;
	}

	std::vector<Intersection> _points{std::begin(points), std::end(points)};

	for (int i = 0; i < _points.size()-1; ++i)
	{
		for (int j = i + 1; j < _points.size(); ++j)
		{
//			Vector3f vBase = _points[j].point - _points[i].point;
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

//				Vector3f vLast = _points[k].point - _points[i].point;
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

	points = std::list<Intersection>{std::begin(_points), std::end(_points)};
}

void VoronoiLattice::OrderPoints(const Vector3f &planeNormal,
								 std::list<Intersection> &points) const
{
	if (!points.empty())
	{
		std::vector<Intersection> _points{std::begin(points), std::end(points)};
//		Point3f pBase = _points.back().point;
		Point3f pBase = _points.back();

		for (int iNext = 0; iNext < _points.size(); ++iNext)
		{
//			Vector3f vBase = _points[iNext].point - pBase;
			Vector3f vBase = _points[iNext] - pBase;

			for (int i = iNext + 1; i <= _points.size(); ++i)
			{
				int iLast = (i == _points.size()) ? 0 : i;
//				Vector3f vLast = _points[iLast].point - pBase;
				Vector3f vLast = _points[iLast] - pBase;
				Vector3f cp = CrossProduct(vBase, vLast);

				if (DotProduct(planeNormal, cp) < EPS_SAME_LINE)
				{
					std::swap(_points[iNext], _points[iLast]);
//					vBase = _points[iNext].point - pBase;
					vBase = _points[iNext] - pBase;
				}
			}

//			pBase = _points[iNext].point;
			pBase = _points[iNext];
		}

		points = std::list<Intersection>{std::begin(_points), std::end(_points)};
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

			Point3f v0ToEdge = CrossProduct(pCurr.normal, edgeV);
//			Point3f v1ToEdge = CrossProduct(edgeV, pBase.normal);

//			bool isOk = false;
			bool isOk = true;
			Vector3f n = pBase.normal;
			n.coordinates[3] = -DotProduct(pBase.center, pBase.normal);
			Point3f x = ProjectPointToPlane(pCurr.center, v0ToEdge, n);
//			Point3f x = IntersectVectors(pCurr.center, pBase.center,
//										 v0ToEdge, v1ToEdge, edgeV, isOk);
			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
				edgeLines.push_back(EdgeLine{edgeV, x});
			}
		}
	}
}

void VoronoiLattice::DefineOrthogonalPlanes(const Cell &baseCell,
											Lattice &lattice,
											std::vector<OrthoPlane> &sitePlanes)
{
	int maxSiteDistance = 2.0*sqrt(2)*lattice.cellSize;
	const auto &ids = baseCell.checkedSiteIds;

	for (int i = 0; i < lattice.sideSize; ++i)
	{
		for (int j = 0; j < lattice.sideSize; ++j)
		{
			for (int k = 0; k < lattice.sideSize; ++k)
			{
				Cell &cell = lattice.cells[i][j][k];

				if (cell.id != baseCell.id)
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
					}
				}
			}
		}
	}

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