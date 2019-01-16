#include "Particle.h"

#include "common.h"
#include "ArgPP.h"
#include "Converter.h"

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

	for (Facet &facet : facets)
	{
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
		ofile << facet /*<< facet.arr[0]*/ << std::endl;
	}

//	center = center/facets.size();
//	ofile << center.point[0] << ' ' << center.point[1] << ' ' << center.point[2] << endl;
	ofile.close();
}

void Triangulate(const std::vector<Facet> &crystal,
				 std::vector<Facet> &triangles)
{
	for (const Facet &facet : crystal)
	{
		if (facet.nVertices == 3) // facet is already triangle
		{
			triangles.push_back(facet);
		}
		else // divide facet into triangles
		{
			for (int i = 1; i+1 < facet.nVertices; ++i)
			{
				Facet triangle;
				triangle.AddVertex(facet.arr[0]); // base vertex
				triangle.AddVertex(facet.arr[i]);
				triangle.AddVertex(facet.arr[i+1]);
				triangles.push_back(triangle);
			}
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

	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
			converter.ReadStl("data/" + filename, triangles);

			if (triangles.empty())
			{
				continue;
			}

			filename = CutSubstring(filename, ".stl");
			filename = dir + filename;
			//	OutputFacets(triangles);

			std::vector<Facet> crystal;
			converter.MergeCrystal(triangles, crystal);

			for (Facet &facet : crystal)
			{
				for (int i = 0; i < facet.nVertices; ++i)
				{
					facet.arr[i] = facet.arr[i]*SIZE_INDEX;
				}
			}

			WriteCry(crystal, filename + "_mbs");
			converter.WriteNat(crystal, filename + "_nat");

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

	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
	return 0;
}
