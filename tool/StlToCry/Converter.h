#pragma once

#include <iostream>
#include <vector>
#include <list>

#include "Facet.h"

class Converter
{
public:
	Converter();
	void ReadStl(const std::string &filename, std::vector<Facet> &triangles);
	void ReadCry(const std::string &filename, std::vector<Facet> &crystal);
	void MergeCrystal(std::vector<Facet> triangles,
					  std::vector<Facet> &mergedFacets);

	void WriteNat(const std::vector<Facet> &crystal,
				  const std::string &outFile);

	void WriteStl(const std::vector<Facet> &triangles,
				  const std::string &outFile);

	void Triangulate(const std::vector<Facet> &crystal,
					 std::vector<Facet> &triangles);
<<<<<<< HEAD
	void Triangulate(const std::vector<std::vector<Point3f>> &facets,
=======
	void Triangulate(const std::vector<std::list<Point3f>> &facets,
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
					 std::vector<Facet> &triangles);
private:
	Point3f ReadVertex(char *buff, char *ptr, char *trash);
	void MergeTriangles(std::vector<Facet> &rest, Facet &convex);
<<<<<<< HEAD
=======
	void FindSameFacetTriangles(std::vector<Facet> &triangles,
								std::vector<Facet> &oneFacetTriangles);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
};
