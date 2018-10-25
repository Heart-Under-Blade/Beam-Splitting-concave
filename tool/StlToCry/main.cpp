#include <iostream>
#include <fstream>
#include <vector>
#include <cstring>
#include <float.h>
#include <limits.h>

#include "Facet.h"
#include "Particle.h"

#define NRM_EPS 10*FLT_EPSILON
#define PNT_EPS FLT_EPSILON

using namespace std;

Point3f ReadVertex(char *buff, char *ptr, char *trash)
{
	Point3f p;
	ptr = strtok(buff, " "); // skip "vertex" word

	ptr = strtok(NULL, " ");
	p.cx = strtod(ptr, &trash);
	ptr = strtok(NULL, " ");
	p.cy = strtod(ptr, &trash);
	ptr = strtok(NULL, " ");
	p.cz = strtod(ptr, &trash);

	return p;
}

void ReadStl(const std::string &filename, std::vector<Facet> &triangles)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	char *ptr, *trash;

	pfile.getline(buff, bufSize); // skip first line "solid HOLDER"
	pfile.getline(buff, bufSize);

	while (std::string(buff) != "endsolid")
	{
		Facet facet;

		// read normal
		ptr = strtok(buff, " "); // skip "facet" word
		ptr = strtok(NULL, " "); // skip "normal" word

		ptr = strtok(NULL, " ");
		facet.ex_normal.cx = strtod(ptr, &trash);
		ptr = strtok(NULL, " ");
		facet.ex_normal.cy = strtod(ptr, &trash);
		ptr = strtok(NULL, " ");
		facet.ex_normal.cz = strtod(ptr, &trash);

		// read vertices
		pfile.getline(buff, bufSize); // skip "outer loop" line

		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));
		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));
		pfile.getline(buff, bufSize);
		facet.AddVertex(ReadVertex(buff, ptr, trash));

		triangles.push_back(facet);

		pfile.getline(buff, bufSize); // skip "endloop" line
		pfile.getline(buff, bufSize); // skip "endfacet" line
		pfile.getline(buff, bufSize);
	}

	pfile.close();
}

int FindEqualPoint(const Point3f &p, const Facet &merged)
{
	if (p.IsEqualTo(merged.arr[0], PNT_EPS))
	{
		return 0;
	}
	else if (p.IsEqualTo(merged.arr[1], PNT_EPS))
	{
		return 1;
	}
	else if (p.IsEqualTo(merged.arr[2], PNT_EPS))
	{
		return 2;
	}
	else // no equal points
	{
		return -1;
	}
}

bool IsNear(int i1, int i2, int n)
{
	return (abs(i2 - i1) == 1) ||
			(i2 == n-1 && i1 == 0) ||
			(i2 == 0 && i1 == n-1);
}

bool FindEqualPoints(const Facet &origin, const Facet &next, Array<int> &points)
{
	int nextIndex;

	for (int i = 0; i < origin.nVertices && points.nElems != 6; ++i)
	{
		nextIndex = FindEqualPoint(origin.arr[i], next);

		if (nextIndex != -1)
		{
			points.Add(i);
			points.Add(nextIndex);
		}
	}

	return (points.nElems >= 4) &&
			IsNear(points.elems[0], points.elems[2], origin.nVertices);
}

void Merge(const Array<int> &points, const Facet &checking, Facet &merged)
{
	if (points.nElems == 4)
	{
		int pointToInsert;

		if ((points.elems[1] == 0 && points.elems[3] == 1) ||
				(points.elems[1] == 1 && points.elems[3] == 0))
		{
			pointToInsert = 2;
		}
		else if ((points.elems[1] == 1 && points.elems[3] == 2) ||
				 (points.elems[1] == 2 && points.elems[3] == 1))
		{
			pointToInsert = 0;
		}
		else //
		{
			pointToInsert = 1;
		}

		int place = (points.elems[0] == merged.nVertices-1) ? 0
															: points.elems[0] + 1;
		merged.InsertVertex(place, checking.arr[pointToInsert]);
#ifdef _DEBUG // DEB
		if (merged.nVertices > 200)
			int gg = 0;
#endif
	}
	else if (points.nElems == 6)
	{
		if (IsNear(points.elems[0], points.elems[2], merged.nVertices) &&
				IsNear(points.elems[2], points.elems[4], merged.nVertices))
		{
			merged.DeleteVertex(points.elems[2]);
		}
	}
	else
	{
		throw std::exception();
	}
}

void OutputFacets(const std::vector<Facet> &facets)
{
	std::string outFile = "out.dat";
	std::ofstream ofile(outFile, std::ios::out);

	if (!ofile.is_open())
	{
		std::cerr << "File \"" << outFile << "\" is not found" << std::endl;
		throw std::exception();
	}

	for (const Facet &facet : facets)
	{
		ofile << facet << std::endl;
	}

	ofile.close();
}

bool IsConvex(const Facet &merged)
{
	double dp1, dp2;
	int i1, i2;
	Vector3f v1 = merged.arr[1] - merged.arr[0];
	Vector3f v2 = merged.arr[2] - merged.arr[0];
	dp1 = Vector3f::DotProduct(v1, v2);

	for (int i = 1; i < merged.nVertices; ++i)
	{
		if (i == merged.nVertices-2)
		{
			i1 = merged.nVertices - 1;
			i2 = 0;
		}
		else if (i == merged.nVertices-1)
		{
			i1 = 0;
			i2 = 1;
		}
		else
		{
			i1 = i + 1;
			i2 = i + 2;
		}

		Vector3f v1 = merged.arr[i1] - merged.arr[i];
		Vector3f v2 = merged.arr[i2] - merged.arr[i];
		dp2 = Vector3f::DotProduct(v1, v2);

		if (dp2*dp1 < 0)
		{
			return false;
		}
	}

	return true;
}

void MergeTriangles(std::vector<Facet> &rest, Facet &convex)
{
	std::vector<Facet> triangles = rest;
	Facet merged = triangles[0];
#ifdef _DEBUG // DEB
	if (triangles.size() == 88)
		int fff = 0;
#endif
	auto it = triangles.begin();
	std::advance(it, 0);
	triangles.erase(it);

	convex = merged;
	rest = triangles;

	int f = 0;
	while (!triangles.empty())
	{
		bool isAdded = false;

		for (int i = 0; i < triangles.size(); ++i)
		{
			const Facet &triangle = triangles[i];
			Array<int> points;

//			if (merged.nVertices == 82)
//			{
//				std::vector<Facet> dd;
//				dd.push_back(merged);
//				OutputFacets(dd);
//				int fff = 0;
//			}

			if (FindEqualPoints(merged, triangle, points))
			{
				Merge(points, triangle, merged);

				it = triangles.begin();
				std::advance(it, i);
				triangles.erase(it);

				if (IsConvex(merged))
				{
					convex = merged;
					rest = triangles;
				}

				isAdded = true;
				break;
			}
		}

		if (!isAdded)
		{
			break;
		}
	}
}

int main()
{
	std::string filename = "particle.stl";
	std::vector<Facet> triangles;
	ReadStl(filename, triangles);

//	OutputFacets(triangles);

	std::vector<Facet> mergedFacets;

	while (!triangles.empty())
	{
		std::vector<Facet> oneFacetTriangles;
		Point3f normal = triangles[0].ex_normal;
		double d1 = -Point3f::DotProduct(triangles[0].arr[0], -normal);

		for (int i = 0; i < triangles.size(); ++i)
		{
			const Facet &checking = triangles[i];

			if (normal.IsEqualTo(checking.ex_normal, NRM_EPS))
			{
				double d2 = -Point3f::DotProduct(checking.arr[0], -normal);

				if (fabs(d1 - d2) < NRM_EPS)
				{
					oneFacetTriangles.push_back(checking);

					auto it = triangles.begin();
					std::advance(it, i);
					triangles.erase(it);
					--i;
				}
			}
		}

		while (!oneFacetTriangles.empty())
		{
			Facet merged;
			//		OutputFacets(oneFacetTriangles);
			MergeTriangles(oneFacetTriangles, merged);
			mergedFacets.push_back(merged);
		}
	}

	OutputFacets(mergedFacets);
	int fff = 0;
}
