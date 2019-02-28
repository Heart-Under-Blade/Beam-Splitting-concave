#include "VoronoiLattice.h"

#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>

#include "geometry_lib.h"
#include "Converter.h"

#define EPS_OUTLINE 0.8001
#define EPS_DUPLICATE 0.003
#define EPS_REDUNDANT 0.003

void VoronoiLattice::OutputFacets(std::vector<std::list<Point3f>> &facets)
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
	std::ofstream file("facets.dat", std::ios::out);

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

VoronoiLattice::VoronoiLattice(double latticeSize, double splitRatio)
	: m_size(latticeSize)
{
	GenerateSites(latticeSize, splitRatio);

	std::vector<Cell> lattice;
#ifdef _DEBUG // DEB
	std::ofstream dfile("debug.dat", std::ios::out);
#endif
	for (int i = /*0*/13; i < 14/*m_sites.size()*/; ++i)
//	for (int i = 0; i < m_sites.size(); ++i)
	{
		std::cout << "Progress: " << i << "/" << m_sites.size() << std::endl;

		DefineFacetPlanes(i);
		Cell facets;

		Site &site = m_sites[i];
#ifdef _DEBUG // DEB
		double k = 10;
		for (auto &p : site.planes)
		{
			dfile << p.center << std::endl;
//			dfile << p.center + p.normal * k << std::endl << std::endl;
		}
		dfile << std::endl;
#endif
		// facets
		for (int j = 0; j < site.planes.size(); ++j)
		{
			std::vector<PointPair> edgeLines;
			DefineFacetEdgeLines(site, j, edgeLines);

#ifdef _DEBUG // DEB
			std::ofstream ffile("facet.dat", std::ios::out);
			for (auto &pp : edgeLines)
			{
				Point3f b = (pp.first * k) + pp.second;
				ffile << b << std::endl;
				ffile << (pp.first * -k) + pp.second << std::endl << std::endl;
			}
//			dfile << site.point << std::endl;
//			dfile << plane.center << std::endl;
			ffile << std::endl;
#endif
			auto &normal = site.planes[j].normal;

			std::list<Point3f> points;
			DefineIntersectionPoints(edgeLines, normal, points);

			RemoveOutlinePoints(i, points);
			RemoveDuplicatedPoints(points);
			RemoveRedundantPoints(points);
			OrderPoints(points, normal);
//			FixPointsOrder(normal, points);

#ifdef _DEBUG // DEB
			ffile << std::endl;
			for (auto &pp : points)
			{
				ffile << pp << std::endl;
			}

			if (!points.empty())
			{
				dfile << site.planes[j].center << std::endl;
				facets.push_back(points);
			}

			ffile.close();
#endif
		}

#ifdef _DEBUG // DEB
		dfile.close();

		OutputFacets(facets);

		Converter converter;
		std::vector<Facet> triangles;
		converter.Triangulate(facets, triangles);
		converter.WriteStl(triangles, "voronoi");
#endif
		lattice.push_back(facets);
	}

//	Converter converter;
//	std::vector<Facet> triangles;

//	for (auto &cell : lattice)
//	{
//		std::vector<Facet> facet;
//		converter.Triangulate(cell, facet);
//		triangles.insert(triangles.end(), facet.begin(), facet.end());
//	}

//	converter.WriteStl(triangles, "voronoi");
}

void VoronoiLattice::DefineIntersectionPoints(
		const std::vector<PointPair> &edgeVectors, const Vector3f &planeNormal,
		std::list<Point3f> &intersectPoints)
{
	for (int i = 0; i < edgeVectors.size(); ++i)
	{
		for (int j = 0; j < edgeVectors.size(); ++j)
		{
			if (j != i)
			{
				bool isOk = false;
				Point3f x = IntersectVectors(
							edgeVectors[i].second, edgeVectors[j].second,
							edgeVectors[i].first, edgeVectors[j].first,
							planeNormal, isOk);
				if (isOk)
				{
					intersectPoints.push_back(x);
				}
			}
		}
	}
}

void VoronoiLattice::RemoveOutlinePoints(int siteIndex,
										 std::list<Point3f> &points)
{
	for (int i = 0; i < m_sites[siteIndex].planes.size(); ++i)
	{
		Point3f &n = m_sites[siteIndex].planes[i].normal;
		Point3f &c = m_sites[siteIndex].planes[i].center;
		points.remove_if([&](Point3f p) {
			return DotProduct(n, p-c) > EPS_OUTLINE;
		});
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
				double len = Length(*it1 - *it0);

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

void VoronoiLattice::RemoveRedundantPoints(std::list<Point3f> &points)
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
				bool isErased = Length(cp) < EPS_REDUNDANT;

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

void VoronoiLattice::OrderPoints(std::list<Point3f> &points,
								 const Vector3f &planeNormal) const
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

				if (DotProduct(planeNormal, cp) < 0)
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

void VoronoiLattice::DefineFacetEdgeLines(const Site &site, int planeIndex,
										  std::vector<PointPair> &edgeLines)
{
	auto &basePlane = site.planes[planeIndex];

	for (int i = 0; i < site.planes.size(); ++i)
	{
		if (i != planeIndex)
		{
			auto &currPlane = site.planes[i];

			Point3f edgeV = CrossProduct(basePlane.normal, currPlane.normal);

			Point3f v0ToEdge = CrossProduct(currPlane.normal, edgeV);
			Point3f v1ToEdge = CrossProduct(edgeV, basePlane.normal);

			bool isOk = false;
			Point3f x = IntersectVectors(currPlane.center, basePlane.center,
										 v0ToEdge, v1ToEdge,
										 edgeV, isOk);
			if (isOk)
			{
				edgeV.coordinates[3] = -DotProduct(x, edgeV);
#ifdef _DEBUG // DEB
				double k = 100;
				Point3f b = (edgeV * k) + x;
				Point3f d = (edgeV * -k) + x;
				Point3f b1 = (v0ToEdge * k) + currPlane.center;
				Point3f d1 = (v1ToEdge * k) + basePlane.center;
				int fff = 0;
#endif
				edgeLines.push_back(PointPair{edgeV, x});
			}
		}
	}
}

void VoronoiLattice::DefineFacetPlanes(int siteIndex)
{
	for (int i = 0; i < m_sites.size(); ++i)
	{
		if (i != siteIndex)
		{
			Point3f normal = m_sites[i].point - m_sites[siteIndex].point;
			Normalize(normal);

			// D-param
			Point3f center = (m_sites[siteIndex].point + m_sites[i].point)/2;
			normal.coordinates[3] = -DotProduct(center, normal);
			m_sites[siteIndex].planes.push_back(OrtoPlane{normal, center});
		}
	}
}

float RandomValueF(float from, float to)
{
	return from + static_cast<float>(rand())/
			(static_cast<float>(RAND_MAX/(to-from)));
}

Point3f RandomPoint(const Point3f &from, const Point3f &to)
{
	Point3f p;
	p.coordinates[0] = RandomValueF(from.coordinates[0], to.coordinates[0]);
	p.coordinates[1] = RandomValueF(from.coordinates[1], to.coordinates[1]);
	p.coordinates[2] = RandomValueF(from.coordinates[2], to.coordinates[2]);
	return p;
}

void VoronoiLattice::GenerateSites(double latticeSize, double splitRatio)
{
	m_sites.resize(splitRatio*splitRatio*splitRatio);
	double cellSize = latticeSize/splitRatio;

	double min = -latticeSize/2;
	double max = latticeSize/2;

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

	int i = 0;
#ifdef _DEBUG // DEB
	std::ofstream file("sites.dat", std::ios::out);
#endif
	for (fromX = min, toX = min + cellSize;
		 toX <= max;
		 fromX += cellSize, toX += cellSize)
	{
		for (fromY = min, toY = min + cellSize;
			 toY <= max;
			 fromY += cellSize, toY += cellSize)
		{
			for (fromZ = min, toZ = min + cellSize;
				 toZ <= max;
				 fromZ += cellSize, toZ += cellSize)
			{
				Point3f p = RandomPoint(from, to);
				m_sites[i++].point = p;
#ifdef _DEBUG // DEB
				file << p.coordinates[0] << " " << p.coordinates[1] << " "
								   << p.coordinates[2] << std::endl;
#endif
			}
		}
	}

#ifdef _DEBUG // DEB
	file.close();
#endif
}
