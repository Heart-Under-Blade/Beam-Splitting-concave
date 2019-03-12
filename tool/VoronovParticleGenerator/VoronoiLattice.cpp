#include "VoronoiLattice.h"

#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>
#include <algorithm>

#include "geometry_lib.h"
#include "Converter.h"

//#define EPS_SAME_FACET 0.003
#define EPS_SAME_FACET 0.43921
#define EPS_DUPLICATE 0.003
#define EPS_SAME_LINE 0.00003

VoronoiLattice::VoronoiLattice(double latticeSize, int splitRatio)
	: m_size(latticeSize)
{
	std::vector<Cell> lattice;
	SiteLattice sites(latticeSize, splitRatio);

#ifdef _DEBUG // DEB
	std::ofstream dfile("debug.dat", std::ios::out);
#endif

	int count = 0;

	for (int i = 1; i < sites.sideSize-1; ++i)
	{
		for (int j = 1; j < sites.sideSize-1; ++j)
		{
			for (int k = 1; k < sites.sideSize-1; ++k)
			{
				std::cout << "Progress: " << count++
						  << "/" << sites.size << std::endl;

				std::vector<OrthoPlane> sitePlanes;
				DefineOrthogonalPlanes(SiteIndex{i, j, k}, sites, sitePlanes);
				Cell cell; // lattice cell consists of facets

#ifdef _DEBUG // DEB
				double kk = 10;
				for (auto &p : sitePlanes)
				{
					dfile << p.center << std::endl;
		//			dfile << p.center + p.normal * kk << std::endl << std::endl;
				}
				dfile << std::endl;
#endif
				for (int planeNo = 0; planeNo < sitePlanes.size(); ++planeNo)
				{
#ifdef _DEBUG // DEB
					if (cell.size() == 8)
					{
						int ass = 0;
					}
#endif
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

					if (cell.size() == 8)
					{
						int ass = 0;

						for (Intersection &is : points)
						{
							if (is.firstPlaneId == 1 || is.secondPlaneId == 1)
								ffile << is.point << std::endl;
						}
					}
#endif
					RemoveExternalPoints(sitePlanes, points);
					RemoveDuplicatedPoints(points);
					RemoveSameLinePoints(points);
					OrderPoints(normal, points);
		//			FixPointsOrder(normal, points);

		//			RemoveDistantPlanes(points, sitePlanes[planeId], sitePlanes);

#ifdef _DEBUG // DEB
					if (!points.empty())
					{
						dfile << sitePlanes[planeNo].center << std::endl;

						std::list<Point3f> facet;

						for (auto &p : points)
						{
							facet.push_back(p.point);
						}

						cell.push_back(facet);
					}
					else
					{
		//				RemovePlane(planeNo, sitePlanes);
		//				--planeNo;
					}

					ffile.close();
#endif
				}

#ifdef _DEBUG // DEB
				dfile.close();
				std::string filename = "cell_" + std::to_string(count);
				OutputFacets(filename + ".dat", cell);

				Converter converter;
				std::vector<Facet> triangles;
				converter.Triangulate(cell, triangles);
				converter.WriteStl(triangles, filename);
#endif
				lattice.push_back(cell);
			}
		}
	}

//	for (auto &cell : lattice)
//	{
//		OutputFacets("cell_" + std::to_string(count), cell);
//	}

	Converter converter;
	std::vector<Facet> triangles;

	for (auto &cell : lattice)
	{
		std::vector<Facet> facet;
		converter.Triangulate(cell, facet);
		triangles.insert(triangles.end(), facet.begin(), facet.end());
	}

	converter.WriteStl(triangles, "voronoi");
}

void VoronoiLattice::RemovePlane(int planeNo, std::vector<OrthoPlane> &sitePlanes)
{
	auto it = sitePlanes.begin();
	std::advance(it, planeNo);
	sitePlanes.erase(it);
}

void VoronoiLattice::OutputFacets(const std::string &filename,
								  std::vector<std::list<Point3f>> &facets)
{
	// find max size of lists
	int maxSize = 0;
	for (int i = 0; i < facets.size(); ++i)
	{
		int newSize = facets[i].size();

		if (newSize > maxSize)
		{
			maxSize = newSize;
		}
	}

	std::vector<std::list<Point3f>::iterator> iters;
	for (int i = 0; i < facets.size(); ++i)
	{
		iters.push_back(facets[i].begin());
	}
	std::ofstream file(filename, std::ios::out);

	for (int i = 0; i < maxSize; ++i)
	{
		for (int j = 0; j < facets.size(); ++j)
		{
			auto &it = iters[j];

			if (it == facets[j].end())
			{
				file << "# # # # ";
			}
			else
			{
				Point3f &p = *it;
				file << p << " # ";
				++it;
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
					if (count == p.firstPlaneId || count == p.secondPlaneId)
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
		for (int j = 0; j < edgeLines.size(); ++j)
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
					points.push_back(Intersection{x, edgeLines[i].planeId,
												  edgeLines[j].planeId});
				}
			}
		}
	}
}

void VoronoiLattice::RemoveExternalPoints(
		const std::vector<OrthoPlane> &sitePlanes,
		std::list<Intersection> &points)
{
	for (int i = 0; i < sitePlanes.size(); ++i)
	{
		const Point3f &n = sitePlanes[i].normal;
		const Point3f &c = sitePlanes[i].center;
		points.remove_if([&](Intersection &p) {
			Vector3f v = p.point-c;
//			Normalize(v);
			double vv = DotProduct(n, v);
			bool res = vv > EPS_SAME_FACET;
			if (res && p.point.coordinates[0] < -40 && p.point.coordinates[0] > -41)
			{
				int fff = 0;
			}
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
				double len = Length((*it1).point - (*it0).point);

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
			Vector3f vBase = _points[j].point - _points[i].point;
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

				Vector3f vLast = _points[k].point - _points[i].point;
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
		Point3f pBase = _points.back().point;

		for (int iNext = 0; iNext < _points.size(); ++iNext)
		{
			Vector3f vBase = _points[iNext].point - pBase;

			for (int i = iNext + 1; i <= _points.size(); ++i)
			{
				int iLast = (i == _points.size()) ? 0 : i;
				Vector3f vLast = _points[iLast].point - pBase;
				Vector3f cp = CrossProduct(vBase, vLast);

				if (DotProduct(planeNormal, cp) < 0)
				{
					std::swap(_points[iNext], _points[iLast]);
					vBase = _points[iNext].point - pBase;
				}
			}

			pBase = _points[iNext].point;
		}

		points = std::list<Intersection>{std::begin(_points), std::end(_points)};
	}
}

void VoronoiLattice::DefineFacetEdgeLines(
		const std::vector<OrthoPlane> &sitePlanes, int planeIndex,
		std::vector<EdgeLine> &edgeLines)
{
	auto &pBase = sitePlanes[planeIndex];

	for (int i = 0; i < sitePlanes.size(); ++i)
	{
		if (i != planeIndex)
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
				edgeLines.push_back(EdgeLine{pCurr.id, edgeV, x});
			}
#ifdef _DEBUG // DEB
			else
			{
				int ff = 0;
			}
#endif
		}
	}
}

void VoronoiLattice::DefineOrthogonalPlanes(const SiteIndex &siteIndex,
											const SiteLattice &sites,
											std::vector<OrthoPlane> &sitePlanes)
{
	for (int i = 0; i < sites.sideSize; ++i)
	{
		for (int j = 0; j < sites.sideSize; ++j)
		{
			for (int k = 0; k < sites.sideSize; ++k)
			{
				if (siteIndex.i != i && siteIndex.j != j && siteIndex.k != k)
				{
					const Point3f &baseSite = sites.GetSite(siteIndex);
					Point3f normal = sites.sites[i][j][k] - baseSite;
					Normalize(normal);

					// D-param
					Point3f center = (baseSite + sites.sites[i][j][k])/2;
					normal.coordinates[3] = -DotProduct(center, normal);
					sitePlanes.push_back(OrthoPlane{i, normal, center});
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
