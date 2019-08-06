#include "Particle.h"

#include "common.h"
#include "ArgPP.h"
#include "Converter.h"

#include <limits.h>
#include <fstream>

#define SIZE_INDEX 1

struct ParticleProperties
{
	double minLen;
	double maxLen;
	double minArea;
	double maxArea;
	int minAreaFacetNo;
	int minLenFacetNo;
	int maxLenFacetNo;
	int maxAreaFacetNo;

	ParticleProperties()
	{
		minLen = LONG_LONG_MAX;
		maxLen = 0;
		minArea = LONG_LONG_MAX;
		maxArea = 0;
		minAreaFacetNo = INT_MAX;
		minLenFacetNo = INT_MAX;
		maxAreaFacetNo = INT_MAX;
		maxLenFacetNo = INT_MAX;
	}
};

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

	ofile << facets.at(0);

	for (int i = 1; i < facets.size(); ++i)
	{
		ofile << std::endl;
		ofile << facets.at(i) /*<< facet.arr[0]*/;
	}

	ofile.close();
}

void AnalyseParticle(const std::vector<Facet> &triangles,
					 ParticleProperties &props)
{
	double newLen;

	for (const Facet &tr : triangles)
	{
		for (int i = 1; i <= tr.nVertices; ++i)
		{
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
				props.maxLenFacetNo = tr.index;
				props.maxLen = newLen;
			}
		}

		double newArea = tr.Area();

		if (newArea < props.minArea)
		{
			props.minAreaFacetNo = tr.index;
			props.minArea = newArea;
		}

		if (newArea > props.maxArea)
		{
			props.maxAreaFacetNo = tr.index;
			props.maxArea = newArea;
		}
	}
}

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
		  << " : min length = " << props.minLen << " (facet " << props.minLenFacetNo
		  << "); max length = " << props.maxLen << " (facet " << props.maxLenFacetNo
		  << "); min area = " << props.minArea << " (facet " << props.minAreaFacetNo
		  << "); max area = " << props.maxArea << " (facet " << props.maxAreaFacetNo
		  << "); Dmax = " << particle->MaximalDimension()
		  << "; particle area = " << particle->Area()
		  << "; volume = " << v
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
	double newArea = -1;

	Converter converter;
	ArgPP parcer;

	if (argc > 1)
	{
		parcer.AddRule("i", 1); // input data type ("cry", "stl")
//		parcer.AddRule("o", '+'); // input data type ("cry", "stl")
		parcer.AddRule("f", '+'); // input files of mask
		parcer.AddRule("rs", 2, true);
//		parcer.AddRule("o", '*'); // output files ("cry", "stl", "nat")

		parcer.Parse(argc, argv);

		auto type = parcer.GetStringValue("i", 0);

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
			else if (type == "a")
			{
				newArea = parcer.GetDoubleValue("rs", 1);
			}
			else
			{
				std::cerr << "Unknown argument \"" << type << "\"" << std::endl;
				throw std::exception();
			}
		}
	}

//	std::string path = ExtractPath(mask);

//	std::vector<std::string> filelist = FindFiles(mask);
	std::vector<std::string> filelist;

	int n = parcer.GetArgNumber("f");

	for (int i = 0; i < n; ++i)
	{
		auto filestr = parcer.GetStringValue("f", i);
		filelist.push_back(filestr);
	}

	std::string	dir = CreateDir("out");

	std::ofstream summaryFile(dir + "summary.txt", std::ios::out);

	if (!summaryFile.is_open())
	{
		std::cerr << "File \"" << "summary.txt" << "\" is not found"
				  << std::endl;
		throw std::exception();
	}

	if (isStl)
	{
		for (auto filename : filelist)
		{
			//std::string filename = "particle.stl";
			std::vector<Facet> triangles;
			converter.ReadStl(filename, triangles);

			filename = CutSubstring(filename, ".stl");
			auto const pos = filename.find_last_of('/');
			const auto file = filename.substr(pos+1);
			filename = dir + file;
			//	OutputFacets(triangles);

			std::vector<Facet> mergedCrystal;
			converter.MergeCrystal(triangles, mergedCrystal);

			// number facets
			int facetNo = 0;

			for (Facet &facet : mergedCrystal)
			{
				facet.index = facetNo++;
			}

//			WriteCry(mergedCrystal, filename + "_mbs");

			Particle *particle = new Particle;
//			particle->SetFromFile(filename + "_mbs.dat"/*, 1,
//					  index*maxWavelength*/);
			particle->SetFromFacets(false, false, mergedCrystal);

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
			else if (newArea > 0)
			{
				double oldA = particle->Area();
				double ratio = newArea/oldA;
				ratio = pow(ratio, 1.0/2);
				particle->Scale(ratio);
			}

			particle->Output(filename + "_mbs.dat");
			ParticleToCrystal(particle, mergedCrystal);

			converter.WriteNat(mergedCrystal, file + "_nat");

			OutputSummary(&summaryFile, mergedCrystal, filename, particle);

			if (triangles.empty())
			{
				continue;
			}

			triangles.clear();
			converter.Triangulate(mergedCrystal, triangles);
			converter.WriteStl(triangles, filename);
		}
	}
	else
	{
		for (auto filename : filelist)
		{
			std::vector<Facet> mergedCrystal;
			std::vector<Facet> triangles;

			Particle *particle = new Particle;
			particle->SetFromFile(filename);

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
			else if (newArea > 0)
			{
				double oldA = particle->Area();
				double ratio = newArea/oldA;
				ratio = pow(ratio, 1.0/2);
				particle->Scale(ratio);
			}

			ParticleToCrystal(particle, mergedCrystal);

			// numer facets
			int facetNo = 0;

			for (Facet &facet : mergedCrystal)
			{
				facet.index = facetNo++;
			}

			OutputSummary(&summaryFile, mergedCrystal, filename, particle);

			converter.Triangulate(mergedCrystal, triangles);
			converter.WriteStl(triangles, filename);
		}
	}

	summaryFile.close();

	std::cout << "Done." << std::endl << "Press <Enter> to exit...";
	getchar();
	return 0;
}
