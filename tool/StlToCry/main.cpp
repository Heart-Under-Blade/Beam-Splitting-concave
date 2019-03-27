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
<<<<<<< HEAD
		  << 180 << ' ' << 360 << std::endl << std::endl;
=======
          << 180 << ' ' << 360 << std::endl << std::endl;
>>>>>>> a35fd73175c41f758864fc5b7ba0285a19c70315

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

<<<<<<< HEAD
=======
bool IsConvex(const Facet &merged)
{
	double dp1, dp2;
	int i1, i2;
	Vector3f v1, v2;

	Vector3f n = merged.Normal();

	v1 = merged.arr[1] - merged.arr[0];
	v2 = merged.arr[2] - merged.arr[0];
	dp1 = Vector3f::DotProduct(Vector3f::CrossProduct(v1, v2), n);

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

		v1 = merged.arr[i1] - merged.arr[i];
		v2 = merged.arr[i2] - merged.arr[i];
		dp2 = Vector3f::DotProduct(Vector3f::CrossProduct(v1, v2), n);

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

	while (!triangles.empty())
	{
		bool isAdded = false;

		for (int i = 0; i < triangles.size(); ++i)
		{
			const Facet &triangle = triangles[i];
			Array<int> points;

			if (FindEqualPoints(merged, triangle, points))
			{
				Merge(points, triangle, merged);

				it = triangles.begin();
				std::advance(it, i);
				triangles.erase(it);

				if (IsConvex(merged))
				{
#ifdef _DEBUG // DEB
					if (merged.nVertices == 26)
					{
//						OutputCrystal(std::vector<Facet>{merged});
						int fff = 0;
					}
#endif
					convex = merged;
					rest = triangles;
				}

				isAdded = true;
			}
		}

		if (!isAdded)
		{
			break;
		}
	}
}

void MergeCrystal(std::vector<Facet> triangles,
				 std::vector<Facet> &mergedFacets)
{
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
}

>>>>>>> a35fd73175c41f758864fc5b7ba0285a19c70315
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

<<<<<<< HEAD
	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
=======
	ofile << "\n}\n";

	ofile.close();
}

int main()
{
	std::string filename = "particle.stl";
	std::vector<Facet> triangles;
	ReadStl(filename, triangles);

//	OutputFacets(triangles);

	std::vector<Facet> crystal;
	MergeCrystal(triangles, crystal);

	WriteCry(crystal);
    WriteNat(crystal);

	triangles.clear();
	Triangulate(crystal, triangles);

	WriteStl(triangles);

>>>>>>> a35fd73175c41f758864fc5b7ba0285a19c70315
	return 0;
}
