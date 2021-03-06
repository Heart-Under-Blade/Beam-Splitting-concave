#include "Particle.h"

<<<<<<< HEAD
//#include "common.h"
=======
#include "common.h"
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
#include "ArgPP.h"
#include "Converter.h"

#include <limits.h>
#include <fstream>

#define SIZE_INDEX 1

<<<<<<< HEAD
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

>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
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
//	Point3f center(0, 0, 0);

//	for (Facet &facet : facets)
//	{
//		center = center + facet.Center();
//	}

//	center = center/facets.size();

=======
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
	ofile << facets.at(0);

	for (int i = 1; i < facets.size(); ++i)
	{
		ofile << std::endl;
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
		ofile << facets.at(i) /*<< facet.arr[0]*/;
	}

	ofile.close();
}

void AnalyseParticle(const std::vector<Facet> &triangles,
					 ParticleProperties &props)
{
	double newLen;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	for (const Facet &tr : triangles)
	{
		for (int i = 1; i <= tr.nVertices; ++i)
		{
<<<<<<< HEAD
			newMin = (i == tr.nVertices)
					? Point3f::Length(tr.arr[0] - tr.arr[i-1])
					: Point3f::Length(tr.arr[i] - tr.arr[i-1]);

			if (newMin < minLen)
			{
				minLenNFacet = tr.index;
				minLen = newMin;
=======
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
			}
		}

		double newArea = tr.Area();

<<<<<<< HEAD
		if (newArea < minArea)
		{
			minNFacet = tr.index;
			minArea = newArea;
=======
		if (newArea < props.minArea)
		{
			props.minAreaFacetNo = tr.index;
			props.minArea = newArea;
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
		}
	}
}

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
void ParticleToCrystal(Particle *particle, std::vector<Facet> &mergedCrystal)
{
	mergedCrystal.clear();

	for (int i = 0; i < particle->nElems; ++i)
	{
		mergedCrystal.push_back(particle->elems[i].original);
	}
}

const string ExtractPath(std::string mask)
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

	if (argc > 1)
	{
		parcer.AddRule("t", 1); // input data type ("cry", "stl")
		parcer.AddRule("f", '+'); // input files of mask
		parcer.AddRule("rs", 2, true);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
//		parcer.AddRule("o", '*'); // output files ("cry", "stl", "nat")

		parcer.Parse(argc, argv);

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

		isStl = (type == "stl");
		isCry = (type == "cry");

		if (parcer.IsCatched("rs"))
		{
			auto type = parcer.GetStringValue("rs", 0);

			if (type == "d")
			{
				newDmax = parcer.GetDoubleValue("rs", 1);
			}
			else if (type == "v")
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
		throw std::exception();
	}

	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
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

			converter.WriteNat(mergedCrystal, file + "_nat");

			OutputSummary(&ofile, mergedCrystal, filename, particle);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

			if (triangles.empty())
			{
				continue;
			}

			triangles.clear();
<<<<<<< HEAD
			Triangulate(crystal, triangles);

=======
			converter.Triangulate(mergedCrystal, triangles);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
			converter.WriteStl(triangles, filename);
		}
	}
	else
	{
		for (auto filename : filelist)
		{
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
			converter.WriteStl(triangles, filename);
		}
	}

	ofile.close();

	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
	return 0;
}
