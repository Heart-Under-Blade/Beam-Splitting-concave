#include "VoronoiLattice.h"

#include <iostream>
#include <fstream>
#include <float.h>
#include <iterator>

#include "geometry_lib.h"

#define EPS_OUTLINE 0.0001
#define EPS_DUPLICATE 0.003
#define EPS_REDUNDANT 0.003

VoronoiLattice::VoronoiLattice(double latticeSize, double splitRatio)
	: m_size(latticeSize)
{
	GenerateSites(latticeSize, splitRatio);

	for (int i = 0/*68*/; i < m_sites.size(); ++i)
	{
		DefineFacetPlanes(i);

#ifdef _DEBUG // DEB
		std::vector<std::list<Point3f>> outpoints;
		std::ofstream file("log.dat", std::ios::out);
#endif
		// facets
		for (int j = 0; j < m_sites[i].planeNormals.size(); ++j)
		{
			std::vector<PointPair> edgeVectors;
			DefineFacetEdgeVectors(m_sites[i], j, edgeVectors);

			std::list<Point3f> points;
			DefineIntersectionPoints(edgeVectors, m_sites[i].planeNormals[j],
									 points);

			RemoveOutlinePoints(i, points);
			RemoveDuplicatedPoints(points);
			RemoveRedundantPoints(points);
			OrderPoints(points, m_sites[i].planeNormals[j]);

#ifdef _DEBUG // DEB
//			for (auto &i : edgeVectors)
//			{
//				file << i.second.point[0] << " " << i.second.point[1] << " "
//								   << i.second.point[2] << std::endl;
//			}
			if (!points.empty())
			{
				outpoints.push_back(points);
//				file << std::endl;
			}
#endif
		}
#ifdef _DEBUG // DEB
		// find max size of lists
		int maxSize = 0;
		for (int i = 0; i < outpoints.size(); ++i)
		{
			int newSize = outpoints[i].size();

			if (newSize > maxSize)
			{
				maxSize = newSize;
			}
		}

		std::vector<std::list<Point3f>::iterator> iters;
		for (int i = 0; i < outpoints.size(); ++i)
		{
			iters.push_back(outpoints[i].begin());
//			bool a = iters[i] == outpoints[i].end();
//			bool b = iters[i] == outpoints[i].begin();
//			int sss = 0;
		}

		// outpur
		for (int i = 0; i < maxSize; ++i)
		{
			for (int j = 0; j < outpoints.size(); ++j)
			{
				auto &it = iters[j];

				if (it == outpoints[j].end())
				{
					file << "# # # # ";
				}
				else
				{
					Point3f &p = *it;
					file << p.point[0] << " "
									   << p.point[1] << " "
									   << p.point[2] << " # ";
					++it;

//					// first point again
//					if (it == outpoints[j].end())
//					{
//						auto it1 = outpoints[j].begin();
//						file << it1->point[0] << " " << it1->point[1] << " "
//											 << it1->point[2] << " # ";
//					}
				}
			}

			file << std::endl;
		}

//		for (auto &i : outpoints)
//		{
//			for (const auto &j : i)
//			{
//				file << j.point[0] << " " << j.point[1] << " "
//								   << j.point[2] << std::endl;
//			}
//		}

		file.close();
#endif
	}
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

				intersectPoints.push_back(x);
			}
		}
	}
}

void VoronoiLattice::RemoveOutlinePoints(int siteIndex,
										 std::list<Point3f> &points)
{
	for (int i = 0; i < m_sites[siteIndex].planeNormals.size(); ++i)
	{
		Point3f &n = m_sites[siteIndex].planeNormals[i];
		Point3f &c = m_sites[siteIndex].centers[i];
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
//	if (!points.empty())
//	{
//		std::vector<Point3f> _points{std::begin(points), std::end(points)};
	Point3f pBase = points.back();

	for (auto itNext = points.begin(); itNext != points.end(); ++itNext)
	{
		Vector3f vBase = *itNext - pBase;

		auto it = itNext;
		++it;

		do
		{
			auto itLast = (it == points.end()) ? points.begin() : it;
			Vector3f vLast = *itLast - pBase;
			Vector3f cp = CrossProduct(vBase, vLast);

			if (Length(cp) < EPS_REDUNDANT)
			{
				if (Length(vBase) < Length(vLast))
				{
					points.erase(itNext);
					vBase = vLast;
				}
				else
				{
					points.erase(itLast);
				}
			}

			++it;
		} while (it != points.end());

		pBase = *itNext;
	}

//		points = std::list<Point3f>{std::begin(_points), std::end(_points)};
//	}

	if (points.size() < 3)
	{
		points.clear();
	}
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

			for (int i = iNext + 1; i < _points.size(); ++i)
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

void VoronoiLattice::DefineFacetEdgeVectors(const Site &site, int planeIndex,
											std::vector<PointPair> &edgeVectors)
{
	for (int i = 0; i < site.planeNormals.size(); ++i)
	{
		if (i != planeIndex)
		{
			Point3f edgeVector = CrossProduct(site.planeNormals[planeIndex],
											  site.planeNormals[i]);
			Point3f vector0ToEdge = CrossProduct(site.planeNormals[i],
												 edgeVector);
			Point3f vector1ToEdge = CrossProduct(edgeVector,
												 site.planeNormals[planeIndex]);
			bool isOk = false;
			Point3f x = IntersectVectors(site.centers[i],
										 site.centers[planeIndex],
										 vector0ToEdge, vector1ToEdge,
										 edgeVector, isOk);
			edgeVector.point[3] = -DotProduct(x, edgeVector);
			edgeVectors.push_back(PointPair{edgeVector, x});
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
			normal.point[3] = -DotProduct(center, normal);

			m_sites[siteIndex].planeNormals.push_back(normal);
			m_sites[siteIndex].centers.push_back(center);
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
	p.point[0] = RandomValueF(from.point[0], to.point[0]);
	p.point[1] = RandomValueF(from.point[1], to.point[1]);
	p.point[2] = RandomValueF(from.point[2], to.point[2]);
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

	float &fromX = from.point[0];
	float &fromY = from.point[1];
	float &fromZ = from.point[2];

	float &toX = to.point[0];
	float &toY = to.point[1];
	float &toZ = to.point[2];

	int i = 0;
#ifdef _DEBUG // DEB
	std::ofstream file("log.dat", std::ios::out);
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
				file << p.point[0] << " " << p.point[1] << " "
								   << p.point[2] << std::endl;
#endif
			}
		}
	}

#ifdef _DEBUG // DEB
	file.close();
#endif
}
