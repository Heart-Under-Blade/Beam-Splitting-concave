#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>
#include <cstring>

#include "CalcTimer.h"
#include "macro.h"
#include "test.h"

#include "Mueller.hpp"

#include "Hexagonal.h"
#include "ConcaveHexagonal.h"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"

#include "global.h"
#include "Beam.h"
#include "PhysMtr.hpp"

#include "Tracer.h"
#include "TracerBackScatterPoint.h"
#include "ArgPP.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

enum class ParticleType : int
{
	Hexagonal = 1,
	Bullet = 2,
	BulletRosette = 3,
	ConcaveHexagonal = 10,
	TiltedHexagonal = 11,
	HexagonalAggregate = 12,
	CertainAggregate = 999
};

Tracks trackGroups;

void ImportTracks(int facetNum)
{
	const int bufSize = 1024;
	ifstream trackFile("tracks2.dat", ios::in);

	if (!trackFile.is_open())
	{
		std::cerr << "Track file not found" << std::endl;
		throw std::exception();
	}

	char *buff = (char*)malloc(sizeof(char) * bufSize);

	TrackGroup buffGroup;

	while (!trackFile.eof())
	{
		trackFile.getline(buff, bufSize);

		vector<int> track;

		char *ptr, *trash;
		ptr = strtok(buff, " ");

		size_t groupIndex = 0;
		bool haveGroup = false;

		while (ptr != NULL)
		{
			if (ptr[0] == ':')
			{
				haveGroup = true;
				ptr = strtok(NULL, " ");
				groupIndex = strtol(ptr, &ptr, 10);

				if (groupIndex >= trackGroups.size())
				{
					for (size_t i = trackGroups.size(); i <= groupIndex; ++i)
					{
						trackGroups.push_back(TrackGroup());
					}
				}

				trackGroups[groupIndex].groupID = groupIndex;
				break;
			}

			int tmp = strtol(ptr, &trash, 10);
			track.push_back(tmp);
			ptr = strtok(NULL, " ");
		}

		BigInteger trackID = 0;

		for (int t : track)
		{
			trackID += (t + 1);
			trackID *= (facetNum + 1);
		}

		if (haveGroup)
		{
			trackGroups[groupIndex].tracks.push_back(track);
			trackGroups[groupIndex].arr[trackGroups[groupIndex].size++] = trackID;
		}
		else
		{
			buffGroup.tracks.push_back(track);
			buffGroup.arr[buffGroup.size++] = trackID;
		}

		track.clear();
	}

	if (buffGroup.size != 0) // добавляем треки без группы в отдельные группы
	{
		for (int i = 0; i < buffGroup.size; ++i)
		{
			TrackGroup newGroup;
			newGroup.arr[newGroup.size++] = buffGroup.arr[i];
			newGroup.tracks.push_back(buffGroup.tracks[i]);
			newGroup.groupID = trackGroups.size();
			trackGroups.push_back(newGroup);
		}
	}
}

void SetArgRules(ArgPP &parser)
{
	int zero = 0;
	parser.AddRule("p", '+'); // particle (type, size, ...)
	parser.AddRule("ri", 1); // reflection index
	parser.AddRule("n", 1); // number of internal reflection
	parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
	parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
	parser.AddRule("go", 0, true); // geometrical optics method
	parser.AddRule("po", 0, true); // phisical optics method
	parser.AddRule("w", 1, true, "po"); // wavelength
	parser.AddRule("b", 2, true, "po"); // beta range (begin, end)
	parser.AddRule("g", 2, true, "po"); // gamma range (begin, end)
	parser.AddRule("conus", 3, true, "po"); // calculate only backscatter cone (radius, phi, theta)
	parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
	parser.AddRule("all", 0, true); // calculate all
	parser.AddRule("o", 1, true); // output file name
	parser.AddRule("close", 0, true); // geometrical optics method
	parser.AddRule("gr", 0, true); // output group files
}

Cone SetCone(ArgPP &parser)
{
	double radius = parser.GetDoubleValue("conus", 0);
	int phiCount = parser.GetDoubleValue("conus", 1);
	int thetaCount = parser.GetDoubleValue("conus", 2);
	return Cone(radius, phiCount, thetaCount);
}

AngleRange GetRange(const ArgPP &parser, const std::string &key,
					Particle *particle)
{
	int number;
	double min, max;

	if (key == "b")
	{
		number = parser.GetIntValue("random", 0);

		if (parser.Catched("b"))
		{
			min = DegToRad(parser.GetDoubleValue(key, 0));
			max = DegToRad(parser.GetDoubleValue(key, 1));
		}
		else
		{
			min = 0;
			max = particle->GetSymmetry().beta;
		}
	}
	else if (key == "g")
	{
		number = parser.GetIntValue("random", 1);

		if (parser.Catched("g"))
		{
			min = DegToRad(parser.GetDoubleValue(key, 0));
			max = DegToRad(parser.GetDoubleValue(key, 1));
		}
		else
		{
			min = 0;
			max = particle->GetSymmetry().gamma;
		}
	}
	else
	{
		cerr << "Error! " << __FUNCTION__;
		throw std::exception();
	}

	return AngleRange(min, max, number);
}

int main(int argc, const char* argv[])
{
	Particle *particle = nullptr;

	bool isCloseConsole = false;

	if (argc <= 1) // has no command line arguments
	{
		cout << endl << "No arguments. Press any key to exit...";
		getchar();
		return 1;
	}

	ArgPP parser;
	SetArgRules(parser);
	parser.Parse(argc, argv);

	// TODO: AggregateBuilder

	double refrIndex = parser.GetDoubleValue("ri");

	if (parser.GetArgNumber("p") == 1)
	{
		std::string filename = parser.GetStringValue("pf");
		particle = new Particle();
		particle->SetFromFile(filename);
		particle->SetRefractiveIndex(complex(refrIndex));
	}
	else
	{
		ParticleType type = (ParticleType)parser.GetIntValue("p", 0);
		double height = parser.GetDoubleValue("p", 1);
		double diameter = parser.GetDoubleValue("p", 2);

		double sup;
		int num;

		switch (type)
		{
		case ParticleType::Hexagonal:
			particle = new Hexagonal(refrIndex, diameter, height);
			break;
		case ParticleType::Bullet:
			sup = (diameter*sqrt(3)*tan(DegToRad(62)))/4;
			particle = new Bullet(refrIndex, diameter, height, sup);
			break;
		case ParticleType::BulletRosette:
			sup = (diameter*sqrt(3)*tan(DegToRad(62)))/4;
			particle = new BulletRosette(refrIndex, diameter, height, sup);
			break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
		case ParticleType::ConcaveHexagonal:
			sup = parser.GetDoubleValue("p", 3);
			particle = new ConcaveHexagonal(refrIndex, diameter, height, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = parser.GetIntValue("p", 3);
			particle = new HexagonalAggregate(refrIndex, diameter, height, num);
			break;
		case ParticleType::CertainAggregate:
			sup = parser.GetDoubleValue("p", 3);
			particle = new CertainAggregate(refrIndex, sup);
			break;
		default:
			assert(false && "ERROR! Incorrect type of particle.");
			break;
		}
	}

	particle->Output();

	int reflNum = parser.GetDoubleValue("n");
	std::string dirName = (parser.Catched("o")) ? parser.GetStringValue("o")
												: "M";
	Tracer tracer(particle, reflNum, dirName);

	if (parser.Catched("gr"))
	{
		tracer.SetIsOutputGroups(true);
	}

	if (parser.Catched("po"))
	{
		ImportTracks(particle->facetNum);

		double wave = parser.GetDoubleValue("w");

		if (parser.Catched("fixed"))
		{
			Cone bsCone = SetCone(parser);

			double beta  = parser.GetDoubleValue("fixed", 0);
			double gamma = parser.GetDoubleValue("fixed", 1);
			tracer.TraceFixedPO(beta, gamma, bsCone, trackGroups, wave);
		}
		else if (parser.Catched("random")) // "random"
		{
			AngleRange beta = GetRange(parser, "b", particle);
			AngleRange gamma = GetRange(parser, "g", particle);

			if (parser.Catched("point"))
			{
				TracerBackScatterPoint tracerBSP(particle, reflNum, dirName);
				tracerBSP.SetIsCalcOther(true);
				tracerBSP.Trace(beta, gamma, trackGroups, wave);
			}
			else
			{
				Cone bsCone = SetCone(parser);

				tracer.TraceRandomPO(beta.number, gamma.number, bsCone,
									 trackGroups, wave);
//				tracer.TraceIntervalPO2(betaR, gammaR, bsCone, trackGroups, wave);
			}
		}
		else
		{
			cout << endl << "error";
		}
	}
	else // go
	{
		int betaR = parser.GetIntValue("random", 0);
		int gammaR = parser.GetIntValue("random", 1);

		if (parser.Catched("all"))
		{
			tracer.TraceRandomGO(betaR, gammaR);
		}
		else
		{
			ImportTracks(particle->facetNum);

			tracer.SetIsCalcOther(true);
			tracer.TraceRandomGO(betaR, gammaR, trackGroups);
		}
		//				tracer.TraceSingleOrGO(45, -90, cellNum, trackGroups);
	}

	if (parser.Catched("close"))
	{
		isCloseConsole = true;
	}

	cout << endl << "done";

	if (!isCloseConsole)
	{
		getchar();
	}

	return 0;
}
