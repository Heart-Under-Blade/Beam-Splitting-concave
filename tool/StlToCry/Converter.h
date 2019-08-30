#pragma once

#include <iostream>
#include <vector>
<<<<<<< HEAD
<<<<<<< HEAD
=======
#include <list>
>>>>>>> feature/track_tree
=======
#include <list>
>>>>>>> origin/feature/voronoi

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
<<<<<<< HEAD
<<<<<<< HEAD
private:
	Point3f ReadVertex(char *buff, char *ptr, char *trash);
	void MergeTriangles(std::vector<Facet> &rest, Facet &convex);
=======

	void Triangulate(const std::vector<Facet> &crystal,
					 std::vector<Facet> &triangles);
	void Triangulate(const std::vector<std::list<Point3f>> &facets,
=======

	void Triangulate(const std::vector<Facet> &crystal,
					 std::vector<Facet> &triangles);
	void Triangulate(const std::vector<std::vector<Point3f>> &facets,
>>>>>>> origin/feature/voronoi
					 std::vector<Facet> &triangles);
private:
	Point3f ReadVertex(char *buff, char *ptr, char *trash);
	void MergeTriangles(std::vector<Facet> &rest, Facet &convex);
<<<<<<< HEAD
	void FindSameFacetTriangles(std::vector<Facet> &triangles,
								std::vector<Facet> &oneFacetTriangles);
>>>>>>> feature/track_tree
=======
>>>>>>> origin/feature/voronoi
};
