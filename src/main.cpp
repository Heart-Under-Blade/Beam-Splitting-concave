
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
#include "TracingConvex.h"
#include "TracingConcave.h"
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

struct OrientationRange
{
	int begin;
	int end;
};

matrix back(4,4),	///< Mueller matrix in backward direction
		forw(4,4);	///< Mueller matrix in forward direction

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

		if (parser.Occured("b"))
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

		if (parser.Occured("g"))
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
	Tracing *tracing = nullptr;

	bool isCloseConsole = false;

	if (argc > 1) // has command line arguments
	{
		ArgPP parser;
		SetArgRules(parser);
		parser.Parse(argc, argv);

		ParticleType pt = (ParticleType)parser.GetIntValue("p", 0);
		double h = parser.GetDoubleValue("p", 1);
		double d = parser.GetDoubleValue("p", 2);

		double ri = parser.GetDoubleValue("ri");
		double sup;
		int num;

		// TODO: AggregateBuilder
		switch (pt)
		{
		case ParticleType::Hexagonal:
			particle = new Hexagonal(ri, d, h);
			break;
		case ParticleType::Bullet:
			sup = parser.GetDoubleValue("p", 3);
			particle = new Bullet(ri, d, h, sup);
			break;
		case ParticleType::BulletRosette:
			sup = (d*sqrt(3)*tan(DegToRad(62)))/4;
			particle = new BulletRosette(ri, d, h, sup);
			break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
		case ParticleType::ConcaveHexagonal:
			sup = parser.GetDoubleValue("p", 3);
			particle = new ConcaveHexagonal(ri, d, h, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = parser.GetIntValue("p", 3);
			particle = new HexagonalAggregate(ri, d, h, num);
			break;
		case ParticleType::CertainAggregate:
			sup = parser.GetDoubleValue("p", 3);
			particle = new CertainAggregate(ri, sup);
			break;
		default:
			assert(false && "ERROR! Incorrect type of particle.");
			break;
		}

//		particle->Move(0, 50, 0);
//		particle->Fix();
//		particle->Rotate(RadToDeg(90), 0, 0);
		particle->Output();

		int reflNum = parser.GetDoubleValue("n");

		if (pt == ParticleType::ConcaveHexagonal ||
				pt == ParticleType::HexagonalAggregate ||
				pt == ParticleType::CertainAggregate||
				pt == ParticleType::BulletRosette)
		{
			tracing = new TracingConcave(particle, incidentDir, isOpticalPath,
										 polarizationBasis, reflNum);
		}
		else
		{
			tracing = new TracingConvex(particle, incidentDir, isOpticalPath,
										polarizationBasis, reflNum);
		}

		std::string dirName = (parser.Occured("o")) ? parser.GetStringValue("o")
												   : "M";
		Tracer tracer(tracing, dirName);

		if (parser.Occured("po"))
		{
			ImportTracks(particle->facetNum);

			double wave = parser.GetDoubleValue("w");

			if (parser.Occured("fixed"))
			{
				Cone bsCone = SetCone(parser);

				double beta  = parser.GetDoubleValue("fixed", 0);
				double gamma = parser.GetDoubleValue("fixed", 1);
				tracer.TraceSingleOrPO(beta, gamma, bsCone, trackGroups, wave);
			}
			else if (parser.Occured("random")) // "random"
			{
				AngleRange beta = GetRange(parser, "b", particle);
				AngleRange gamma = GetRange(parser, "g", particle);

				if (parser.Occured("point"))
				{
					tracer.setIsCalcOther(true);
					tracer.TraceBackScatterPointPO(beta, gamma, trackGroups, wave);
				}
				else
				{
					Cone bsCone = SetCone(parser);

					tracer.TraceRandomPO(beta.number, gamma.number, bsCone,
										 trackGroups, wave);
//					tracer.TraceIntervalPO2(betaR, gammaR, bsCone, trackGroups, wave);
				}
			}
			else
			{
				cout << endl << "error";
			}
		}
		else if (parser.Occured("go"))
		{
			int betaR = parser.GetIntValue("random", 0);
			int gammaR = parser.GetIntValue("random", 1);

			if (parser.Occured("all"))
			{
				tracer.TraceIntervalGO(betaR, gammaR);
			}
			else
			{
				ImportTracks(particle->facetNum);

				tracer.setIsCalcOther(true);
				tracer.TraceIntervalGO(betaR, gammaR, trackGroups);
			}
			// tracer.TraceSingleOrGO(45, -90, cellNum, trackGroups);
		}
		else
		{
			cout << endl << "error";
		}

		cout << endl << "done";

		if (parser.Occured("close"))
		{
			isCloseConsole = true;
		}
	}

	if (!isCloseConsole)
	{
		getchar();
	}

	return 0;
}
