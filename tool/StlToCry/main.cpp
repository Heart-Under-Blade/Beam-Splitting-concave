#include "Particle.h"

//#include "common.h"
#include "ArgPP.h"
#include "Converter.h"

#include <limits.h>
#include <fstream>

#define SIZE_INDEX 1

using namespace std;

void WriteCry(std::vector<Facet> &facets, const std::string &outFile)
{
	std::ofstream ofile(outFile + ".dat", std::ios::out);

	if (!ofile.is_open())
	{
		std::cerr << "File \"" << outFile << "\" is not found" << std::endl;
		throw std::exception();
	}

	ofile << 0 << std::endl
		  << 0 << std::endl
		  << 180 << ' ' << 360 << std::endl << std::endl;

//	Point3f center(0, 0, 0);

//	for (Facet &facet : facets)
//	{
//		center = center + facet.Center();
//	}

//	center = center/facets.size();

	ofile << facets.at(0);

	for (int i = 1; i < facets.size(); ++i)
	{
		ofile << std::endl;
//	for (Facet &facet : facets)
//    {
//		const Point3f n = facet.Normal();

//		if (Point3f::DotProduct(facet.Center()-center, n) > 0)
//		{
//			for (int i = 0; i < facet.nVertices/2; ++i)
//			{
//				Point3f buf = facet.arr[i];
//				facet.arr[i] = facet.arr[(facet.nVertices-1)-i];
//				facet.arr[(facet.nVertices-1)-i] = buf;
//			}
//		}

//		ofile << center.point[0] << ' ' << center.point[1] << ' ' << center.point[2] << endl
//								 << n.point[0] << ' '
//								 << n.point[1] << ' '
//								 << n.point[2] << endl ;
#ifdef _DEBUG // DEB
//		if (facet.Area() < 3)
//			continue;
#endif
		ofile << facets.at(i) /*<< facet.arr[0]*/;
	}

//	center = center/facets.size();
//	ofile << center.point[0] << ' ' << center.point[1] << ' ' << center.point[2] << endl;
	ofile.close();
}

void Verify(const std::vector<Facet> &triangles,
			double &minLen, double &minArea, int &minNFacet, int &minLenNFacet)
{
	minLen = LONG_LONG_MAX;
	minArea = LONG_LONG_MAX;
	minNFacet = INT_MAX;
	minLenNFacet = INT_MAX;
	double newMin;

	for (const Facet &tr : triangles)
	{
		for (int i = 1; i <= tr.nVertices; ++i)
		{
			newMin = (i == tr.nVertices)
					? Point3f::Length(tr.arr[0] - tr.arr[i-1])
					: Point3f::Length(tr.arr[i] - tr.arr[i-1]);

			if (newMin < minLen)
			{
				minLenNFacet = tr.index;
				minLen = newMin;
			}
		}

		double newArea = tr.Area();

		if (newArea < minArea)
		{
			minNFacet = tr.index;
			minArea = newArea;
		}
	}
}

int main(int argc, const char *argv[])
{
	bool isStl = true;
	bool isCry = false;
	std::string mask = "data/*.stl";
//	std::string mask = "data/*.dat";
	Converter converter;

	if (argc > 1)
	{
		ArgPP parcer;
		parcer.AddRule("i", 2); // input data type ("cry", "stl") and dir
//		parcer.AddRule("o", '*'); // output files ("cry", "stl", "nat")

		parcer.Parse(argc, argv);

		auto type = parcer.GetStringValue("i", 1);
		mask = parcer.GetStringValue("i", 2);

		isStl = (type == "stl");
		isCry = (type == "cry");
	}

	std::vector<std::string> filelist = FindFiles(mask);
	std::string	dir = CreateDir("out");

	std::ofstream ofile("data/files.dat", std::ios::out);
	if (!ofile.is_open())
	{
		std::cerr << "File \"" << "data/files.dat" << "\" is not found" << std::endl;
		throw std::exception();
	}

	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
			converter.ReadStl("data/" + filename, triangles);

			filename = CutSubstring(filename, ".stl");
			filename = dir + filename;
			//	OutputFacets(triangles);

			std::vector<Facet> crystal;
			converter.MergeCrystal(triangles, crystal);

			int facetNo = 0;

			for (Facet &facet : crystal)
			{
				for (int i = 0; i < facet.nVertices; ++i)
				{
					facet.arr[i] = facet.arr[i]*SIZE_INDEX;
				}

				facet.index = facetNo++;
			}

			WriteCry(crystal, filename + "_mbs");
			converter.WriteNat(crystal, filename + "_nat");

			double minLen;
			double minArea;
			int minLenNFacet;
			int minNFacet;
			Verify(crystal, minLen, minArea, minNFacet, minLenNFacet);

			Particle *particle = new Particle;
			particle->SetFromFile(filename + "_mbs.dat");

			ofile << filename
				  << " : min length = " << minLen
				  << "(facet " << minLenNFacet << ")"
				  << ", min area = " << minArea
				  << "(facet " << minNFacet
				  << "), Dmax = " << particle->ComputeDmax()
				  << std::endl;

			if (triangles.empty())
			{
				continue;
			}

			triangles.clear();
			Triangulate(crystal, triangles);

			converter.WriteStl(triangles, filename);
		}
	}
	else
	{
		for (auto filename : filelist)
		{
			std::vector<Facet> crystal;
			std::vector<Facet> triangles;

			Particle *particle = new Particle;
			particle->SetFromFile("data/" + filename);

			for (int i = 0; i < particle->nElems; ++i)
			{
				crystal.push_back(particle->elems[i].origin);
			}

			Triangulate(crystal, triangles);
			converter.WriteStl(triangles, filename);
		}
	}

	ofile.close();

	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
	return 0;
}
