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

#include "global.h"
#include "Beam.h"
#include "PhysMtr.hpp"

#include "ArgParser.h"
#include "Tracer.h"
#include "TracingConvex.h"
#include "TracingConcave.h"

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
	ConcaveHexagonal = 10,
	TiltedHexagonal = 11,
	HexagonalAggregate = 12,
	CertainAggregate = 999
};

struct OrientationRange
{
	int begin;
	int end;
};

matrix back(4,4),	///< Mueller matrix in backward direction
		forw(4,4);	///< Mueller matrix in forward direction

Arr2D mxd(0, 0, 0, 0);
double sizeBin;
double gammaNorm, betaNorm;
double incomingEnergy;
Point3f incidentDir(0, 0, -1);
Point3f polarizationBasis(0, 1, 0); ///< Basis for polarization characteristic of light
int assertNum = 0;

bool isPhisOptics = true;
bool isOpticalPath = true;
vector<Arr2DC> J; // Jones matrices
double df = 0, dt = 0;
unsigned int MuellerMatrixNumber = 0;
unsigned int NumSum = 16;
int betaMax = 3;

Tracks trackGroups;

int groupCount = 0;

void ImportTracks(int facetNum)
{
	const int bufSize = 1024;
	ifstream trackFile("tracks2.dat", ios::in);
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	TrackGroup buffGroup;

	while (!trackFile.eof())
	{
		trackFile.getline(buff, bufSize);

		vector<int> track;

		char *ptr, *trash;
		ptr = strtok(buff, " ");

		int groupIndex = 0;
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
					for (int i = trackGroups.size();
						 i <= groupIndex; ++i)
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

		long long int trackID = 0;

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

// TODO: написать свой ArgParser

void setAvalableArgs(ArgParser &parser)
{
	parser.addArgument("-p", "--particle", '+', false);
	parser.addArgument("--ri", 1, false);
	parser.addArgument("-n", "--interReflNum", 1, false);
	parser.addArgument("--point", 2); // REF: заменить на --fixed
	parser.addArgument("--range"); // REF: перенести бета и гамма сюда в кач-ве параметров | заменить на --random
	parser.addArgument("-b", "--beta", 3); // REF: попробовать задавать в частице
	parser.addArgument("-g", "--gamma", 3);
	parser.addArgument("-t", "--cellCount", 1);
	parser.addArgument("--po"); // REF: ввести --go
	parser.addArgument("-w", "--wavelength", 1);
	parser.addArgument("--conus", 3);
	parser.addArgument("-o", "--output", 1);
	parser.addArgument("--all");
}

AngleRange GetRange(const char *name, double normCoef, ArgParser &parser)
{
	vector<string> argInterval = parser.retrieve<vector<string>>(name);
	double begin = parser.argToValue<double>(argInterval[0]);
	double end = parser.argToValue<double>(argInterval[1]);
	double count = parser.argToValue<int>(argInterval[2]);
	AngleRange range(begin, end, count, normCoef);
	return range;
}

Cone SetCone(ArgParser &parser)
{
	vector<string> cone_arg = parser.retrieve<vector<string>>("conus");
	double radius = parser.argToValue<double>(cone_arg[0]);
	int phiCount = parser.argToValue<double>(cone_arg[1]);
	int thetaCount = parser.argToValue<double>(cone_arg[2]);
	return Cone(radius, phiCount, thetaCount);
}

int main(int argc, const char** argv)
{
	Particle *particle = nullptr;
	Tracing *tracing = nullptr;

	if (argc > 1) // has command line arguments
	{
		ArgParser parser;
		setAvalableArgs(parser);
		parser.parse(argc, argv);

		vector<string> vec = parser.retrieve<vector<string>>("particle");
		ParticleType pt = (ParticleType)parser.argToValue<int>(vec[0]);
		double h = parser.argToValue<double>(vec[1]);
		double d = parser.argToValue<double>(vec[2]);

		double ri = parser.getArgValue<double>("ri");
		double sup;
		int num;

		// TODO: AggregateBuilder
		switch (pt)
		{
		case ParticleType::Hexagonal:
			particle = new Hexagonal(ri, d, h);
			break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
		case ParticleType::ConcaveHexagonal:
			sup = parser.argToValue<double>(vec[3]);
			particle = new ConcaveHexagonal(ri, d, h, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = parser.argToValue<int>(vec[3]);
			particle = new HexagonalAggregate(ri, d, h, num);
			break;
		case ParticleType::CertainAggregate:
			num = parser.argToValue<int>(vec[3]);
			particle = new CertainAggregate(ri, num);
			break;
		default:
			assert(false && "ERROR! Incorrect type of particle.");
			break;
		}

		toFile(*particle);

		int reflNum = parser.getArgValue<double>("n");

		if (pt == ParticleType::ConcaveHexagonal ||
				pt == ParticleType::HexagonalAggregate)
		{
			tracing = new TracingConcave(particle, incidentDir, isOpticalPath,
										 polarizationBasis, reflNum);
		}
		else
		{
			tracing = new TracingConvex(particle, incidentDir, isOpticalPath,
										polarizationBasis, reflNum);
		}

		double wave = parser.getArgValue<double>("wavelength");

//		ImportTracks(particle->facetNum);
		ImportTracks(particle->facetNum);
		Tracer tracer(tracing, "M");

		if (parser.count("--po") != 0)
		{
			Cone bsCone = SetCone(parser);

			if (parser.count("point") != 0)
			{
				vector<string> vec = parser.retrieve<vector<string>>("point");
				double beta  = parser.argToValue<double>(vec[0]);
				double gamma = parser.argToValue<double>(vec[1]);
				//			beta = 32; gamma = 30;
				tracer.TraceSingleOrPO(beta, gamma, bsCone, trackGroups, wave);
			}
			else // "range"
			{
				AngleRange betaR = GetRange("beta", particle->GetSymmetryBeta(), parser);
				AngleRange gammaR = GetRange("gamma", particle->GetSymmetryGamma(), parser);
				tracer.TraceRandomPO(betaR, gammaR, bsCone, trackGroups, wave);
	//			tracer.TraceIntervalPO2(betaR, gammaR, bsCone, trackGroups, wave);
			}
		}
		else
		{
			if (parser.count("t") != 0)
			{
				AngleRange betaR = GetRange("beta", particle->GetSymmetryBeta(), parser);
				AngleRange gammaR = GetRange("gamma", particle->GetSymmetryGamma(), parser);

				int cellNum = parser.getArgValue<int>("t");

				if (parser.count("all") != 0)
				{
					tracer.TraceIntervalGO(betaR, gammaR, cellNum);
				}
				else
				{
					tracer.setIsCalcOther(true);
//					tracer.TraceBackScatterPointPO(betaR, gammaR, trackGroups, 0.532);
					tracer.TraceIntervalGO(betaR, gammaR, cellNum, trackGroups);
				}
				//				tracer.TraceSingleOrGO(45, -90, cellNum, trackGroups);
			}
			else
			{
				cout << endl << "error";
			}
		}

		cout << endl << "done";
	}

	getchar();
	return 0;
}
