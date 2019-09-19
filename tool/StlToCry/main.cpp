#include "Particle.h"

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
//#include "common.h"
=======
#include "common.h"
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
#include "common.h"
>>>>>>> origin/refactor
=======
#include "common.h"
>>>>>>> origin/refactor
#include "ArgPP.h"
#include "Converter.h"
<<<<<<< HEAD

#include <fstream>

#define SIZE_INDEX 1
=======

#include <limits.h>
#include <fstream>
>>>>>>> feature/track_tree

#define SIZE_INDEX 1

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
=======
void WriteCry(std::vector<Facet> &facets, const std::string &outFile)
>>>>>>> origin/refactor
=======
void WriteCry(std::vector<Facet> &facets, const std::string &outFile)
>>>>>>> origin/refactor
=======
struct ParticleProperties
{
	double minLen;
	double maxLen;
	double minArea;
	int minAreaFacetNo;
	int minLenFacetNo;

	ParticleProperties()
	{
		minLen = LONG_LONG_MAX;
		maxLen = 0;
		minArea = LONG_LONG_MAX;
		minAreaFacetNo = INT_MAX;
		minLenFacetNo = INT_MAX;
	}
};

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
=======
//#include "common.h"
#include "ArgPP.h"
#include "Converter.h"

#include <limits.h>
#include <fstream>

#define SIZE_INDEX 1

>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
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

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
//	Point3f center(0, 0, 0);

//	for (Facet &facet : facets)
//	{
//		center = center + facet.Center();
//	}

//	center = center/facets.size();

=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
	ofile << facets.at(0);

	for (int i = 1; i < facets.size(); ++i)
	{
		ofile << std::endl;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		ofile << facets.at(i) /*<< facet.arr[0]*/;
	}

	ofile.close();
}

void AnalyseParticle(const std::vector<Facet> &triangles,
					 ParticleProperties &props)
{
	double newLen;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor

	for (const Facet &tr : triangles)
	{
		for (int i = 1; i <= tr.nVertices; ++i)
		{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			newMin = (i == tr.nVertices)
					? Point3f::Length(tr.arr[0] - tr.arr[i-1])
					: Point3f::Length(tr.arr[i] - tr.arr[i-1]);

			if (newMin < minLen)
			{
				minLenNFacet = tr.index;
				minLen = newMin;
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			newLen = (i == tr.nVertices)
					? Point3f::Length(tr.vertices[0] - tr.vertices[i-1])
					: Point3f::Length(tr.vertices[i] - tr.vertices[i-1]);

			if (newLen < props.minLen)
			{
				props.minLenFacetNo = tr.index;
				props.minLen = newLen;
			}

			if (newLen > props.maxLen)
			{
				props.maxLen = newLen;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			}
		}

		double newArea = tr.Area();

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
		if (newArea < minArea)
		{
			minNFacet = tr.index;
			minArea = newArea;
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		if (newArea < props.minArea)
		{
			props.minAreaFacetNo = tr.index;
			props.minArea = newArea;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		}
	}
}

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
void ParticleToCrystal(Particle *particle, std::vector<Facet> &mergedCrystal)
>>>>>>> feature/track_tree
{
	mergedCrystal.clear();

	for (int i = 0; i < particle->nElems; ++i)
	{
		mergedCrystal.push_back(particle->elems[i].original);
	}
<<<<<<< HEAD

	ofile << 0 << std::endl
		  << 0 << std::endl
<<<<<<< HEAD
		  << 180 << ' ' << 360 << std::endl << std::endl;
=======
          << 180 << ' ' << 360 << std::endl << std::endl;
>>>>>>> a35fd73175c41f758864fc5b7ba0285a19c70315

=======
>>>>>>> origin/feature/voronoi
//	Point3f center(0, 0, 0);

//	for (Facet &facet : facets)
//	{
//		center = center + facet.Center();
//	}

//	center = center/facets.size();

<<<<<<< HEAD
	for (Facet &facet : facets)
    {
=======
	ofile << facets.at(0);

	for (int i = 1; i < facets.size(); ++i)
	{
		ofile << std::endl;
//	for (Facet &facet : facets)
//    {
>>>>>>> origin/feature/voronoi
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
<<<<<<< HEAD
		ofile << facet /*<< facet.arr[0]*/ << std::endl;
=======
		ofile << facets.at(i) /*<< facet.arr[0]*/;
>>>>>>> origin/feature/voronoi
	}

//	center = center/facets.size();
//	ofile << center.point[0] << ' ' << center.point[1] << ' ' << center.point[2] << endl;
	ofile.close();
}

<<<<<<< HEAD
<<<<<<< HEAD
=======
bool IsConvex(const Facet &merged)
=======
}

const string ExtractPath(std::string mask)
>>>>>>> feature/track_tree
{
	auto pos = mask.find_last_of('/');

	if (pos == std::string::npos)
	{
		pos = mask.find_last_of('\\');

		if (pos == std::string::npos)
		{
			std::cerr << "Path \"" << mask << "\" is not found" << std::endl;
			throw std::exception();
		}
	}

	return mask.substr(0, pos+1);
}

void OutputSummary(std::ofstream *ofile,
				   const std::vector<Facet> &mergedCrystal,
				   const std::string &filename,
				   Particle *particle)
{
	ParticleProperties props;
	AnalyseParticle(mergedCrystal, props);

	double v = particle->Volume();
	double rEq = pow((3*v)/4*M_PI, 1.0/3);

	(*ofile) << filename
		  << " : min length = " << props.minLen
		  << " (facet " << props.minLenFacetNo
		  << "); max length = " << props.maxLen
		  << "; min area = " << props.minArea
		  << " (facet " << props.minAreaFacetNo
		  << "); Dmax = " << particle->MaximalDimension()
		  << "; particle: area = " << particle->Area()
		  << ", volume = " << v
		  << "; R-eq = " << rEq
		  << std::endl;
}

int main(int argc, const char *argv[])
{
	double maxWavelength = 1.064;
	double index = 5;

	bool isStl = true;
	bool isCry = false;
	std::string mask = "data/*.stl";
	double newDmax = -1;
	double newVolume = -1;

	Converter converter;
	ArgPP parcer;

<<<<<<< HEAD
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
=======
	if (argc > 1)
	{
		parcer.AddRule("t", 1); // input data type ("cry", "stl")
		parcer.AddRule("f", '+'); // input files of mask
		parcer.AddRule("rs", 2, true);
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
//		parcer.AddRule("o", '*'); // output files ("cry", "stl", "nat")
>>>>>>> feature/track_tree

		parcer.Parse(argc, argv);

		auto type = parcer.GetStringValue("t", 0);

<<<<<<< HEAD
void MergeCrystal(std::vector<Facet> triangles,
				 std::vector<Facet> &mergedFacets)
{
	while (!triangles.empty())
	{
		std::vector<Facet> oneFacetTriangles;
		Point3f normal = triangles[0].ex_normal;
		double d1 = -Point3f::DotProduct(triangles[0].arr[0], -normal);
=======
		isStl = (type == "stl");
		isCry = (type == "cry");
>>>>>>> feature/track_tree

		if (parcer.IsCatched("rs"))
		{
			auto type = parcer.GetStringValue("rs", 0);

<<<<<<< HEAD
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
=======
			if (type == "d")
			{
				newDmax = parcer.GetDoubleValue("rs", 1);
			}
			else if (type == "v")
>>>>>>> feature/track_tree
			{
				newVolume = parcer.GetDoubleValue("rs", 1);
			}
			else
			{
				std::cerr << "Unknown argument \"" << type << "\"" << std::endl;
				throw std::exception();
			}
		}
	}

<<<<<<< HEAD
=======
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

>>>>>>> origin/feature/voronoi
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
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
//		parcer.AddRule("o", '*'); // output files ("cry", "stl", "nat")

		parcer.Parse(argc, argv);

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
		auto type = parcer.GetStringValue("t", 0);
=======
		auto type = parcer.GetStringValue("i", 1);
		mask = parcer.GetStringValue("i", 2);
>>>>>>> origin/refactor
=======
		auto type = parcer.GetStringValue("i", 1);
		mask = parcer.GetStringValue("i", 2);
>>>>>>> origin/refactor

		isStl = (type == "stl");
		isCry = (type == "cry");
	}

	std::vector<std::string> filelist = FindFiles(mask);
	std::string	dir = CreateDir("out");

<<<<<<< HEAD
=======
	std::ofstream ofile("data/files.dat", std::ios::out);
	if (!ofile.is_open())
	{
		std::cerr << "File \"" << "data/files.dat" << "\" is not found" << std::endl;
		throw std::exception();
	}

>>>>>>> origin/feature/voronoi
	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
			converter.ReadStl("data/" + filename, triangles);

<<<<<<< HEAD
			if (triangles.empty())
			{
				continue;
			}

=======
>>>>>>> origin/feature/voronoi
			filename = CutSubstring(filename, ".stl");
			filename = dir + filename;
			//	OutputFacets(triangles);

			std::vector<Facet> crystal;
			converter.MergeCrystal(triangles, crystal);

<<<<<<< HEAD
=======
			int facetNo = 0;

>>>>>>> origin/feature/voronoi
			for (Facet &facet : crystal)
			{
				for (int i = 0; i < facet.nVertices; ++i)
				{
					facet.arr[i] = facet.arr[i]*SIZE_INDEX;
				}
<<<<<<< HEAD
=======

				facet.index = facetNo++;
>>>>>>> origin/feature/voronoi
			}

			WriteCry(crystal, filename + "_mbs");
			converter.WriteNat(crystal, filename + "_nat");

<<<<<<< HEAD
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
=======
//	std::string path = ExtractPath(mask);
	std::string path = "";

//	std::vector<std::string> filelist = FindFiles(mask);
	std::vector<std::string> filelist;

	int n = parcer.GetArgNumber("f");

	for (int i = 0; i < n; ++i)
	{
		auto filestr = parcer.GetStringValue("f", i);
		filelist.push_back(filestr);
	}

	std::string	dir = CreateDir("out");

	std::ofstream ofile(path + "summary.txt", std::ios::out);

	if (!ofile.is_open())
	{
		std::cerr << "File \"" << path + "summary.txt" << "\" is not found"
				  << std::endl;
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
		throw std::exception();
	}

	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			converter.ReadStl(path + filename, triangles);

			filename = CutSubstring(filename, ".stl");
			auto const pos = filename.find_last_of('/');
			const auto file = filename.substr(pos+1);
			filename = dir + file;
			//	OutputFacets(triangles);

			std::vector<Facet> mergedCrystal;
			converter.MergeCrystal(triangles, mergedCrystal);

			// numer facets
			int facetNo = 0;

			for (Facet &facet : mergedCrystal)
			{
				facet.index = facetNo++;
			}

			WriteCry(mergedCrystal, filename + "_mbs");

			Particle *particle = new Particle;
			particle->SetFromFile(filename + "_mbs.dat"/*, 1,
					  index*maxWavelength*/);

			if (newDmax > 0)
			{
				double oldDMax = particle->MaximalDimension();
				double ratio = newDmax/oldDMax;
				particle->Scale(ratio);
			}
			else if (newVolume > 0)
			{
				double oldV = particle->Volume();
				double ratio = newVolume/oldV;
				ratio = pow(ratio, 1.0/3);
				particle->Scale(ratio);
			}

			ParticleToCrystal(particle, mergedCrystal);
>>>>>>> feature/track_tree

			converter.WriteNat(mergedCrystal, file + "_nat");

<<<<<<< HEAD
<<<<<<< HEAD
			OutputSummary(&ofile, mergedCrystal, filename, particle);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
<<<<<<< HEAD
=======
>>>>>>> origin/refactor
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
=======
			OutputSummary(&ofile, mergedCrystal, filename, particle);
=======
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
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor

			if (triangles.empty())
			{
				continue;
			}

			triangles.clear();
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			Triangulate(crystal, triangles);

=======
			converter.Triangulate(mergedCrystal, triangles);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
			converter.Triangulate(mergedCrystal, triangles);
=======
			Triangulate(crystal, triangles);

>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			converter.WriteStl(triangles, filename);
		}
	}
	else
	{
		for (auto filename : filelist)
		{
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
			std::vector<Facet> crystal;
			std::vector<Facet> triangles;

			Particle *particle = new Particle;
			particle->SetFromFile("data/" + filename);

			for (int i = 0; i < particle->nElems; ++i)
			{
				crystal.push_back(particle->elems[i].origin);
			}

			Triangulate(crystal, triangles);
=======
=======
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			std::vector<Facet> mergedCrystal;
			std::vector<Facet> triangles;

			Particle *particle = new Particle;
			particle->SetFromFile(path + filename);

			if (newDmax > 0)
			{
				double oldDMax = particle->MaximalDimension();
				double ratio = oldDMax/newDmax;
				particle->Scale(ratio);
			}
			else if (newVolume > 0)
			{
				double oldV = particle->Volume();
				double ratio = oldV/newVolume;
				pow(ratio, 1.0/3);
				particle->Scale(ratio);
			}

			ParticleToCrystal(particle, mergedCrystal);

			// numer facets
			int facetNo = 0;

			for (Facet &facet : mergedCrystal)
			{
				facet.index = facetNo++;
			}

			ParticleProperties props;
			AnalyseParticle(mergedCrystal, props);

			OutputSummary(&ofile, mergedCrystal, filename, particle);

			converter.Triangulate(mergedCrystal, triangles);
<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
=======
=======
>>>>>>> origin/refactor
=======
			std::vector<Facet> crystal;
			std::vector<Facet> triangles;

			Particle *particle = new Particle;
			particle->SetFromFile("data/" + filename);

			for (int i = 0; i < particle->nElems; ++i)
			{
				crystal.push_back(particle->elems[i].origin);
			}

			Triangulate(crystal, triangles);
>>>>>>> origin/feature/voronoi
<<<<<<< HEAD
>>>>>>> origin/refactor
=======
>>>>>>> origin/refactor
			converter.WriteStl(triangles, filename);
		}
	}

	ofile.close();

	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
<<<<<<< HEAD
>>>>>>> feature/track_tree
=======
>>>>>>> origin/feature/voronoi
	return 0;
}
