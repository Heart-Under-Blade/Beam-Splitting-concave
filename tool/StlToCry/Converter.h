#pragma once

#include <iostream>
#include <vector>

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
private:
	Point3f ReadVertex(char *buff, char *ptr, char *trash);
	void MergeTriangles(std::vector<Facet> &rest, Facet &convex);
};