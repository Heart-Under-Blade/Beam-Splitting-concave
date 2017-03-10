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

#include "TracingConcave.h"
#include "TracingConvex.h"

#include "global.h"
#include "Beam.h"
#include "PhysMtr.hpp"

#include "argparse.hpp"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

struct OrientationRange
{
	int begin;
	int end;
};

struct Cone
{
	double radius;
	int phi;
	int theta;
};

enum class ParticleType : int
{
	Hexagonal = 1,
	ConcaveHexagonal = 10,
	TiltedHexagonal = 11
};

struct CLArguments
{
	ParticleType particleType;
	double halfHeight;
	double radius;
	double cavityDepth;
	double tiltAngle;
	complex refractionIndex;
	OrientationRange betaRange;
	OrientationRange gammaRange;
	int thetaNumber;
	int interReflNum;
	bool isRandom = false;
	string outfile;
	Cone bsCone; ///< конус в направлении назад
	double wavelength;
} params;

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

struct TrackGroup
{
	int groupID;
	long long int arr[1024];
	int size = 0;
}
trackGroups[32];

int groupCount = 0;

void HandleBeams(vector<Beam> &outBeams, double betaDistrProb, const Tracing &tracer);
void ExtractPeaks(int EDF, double NRM, int ThetaNumber);

void WriteResultsToFile(int ThetaNumber, double NRM, const string &filename);
void WriteStatisticsToConsole(int orNumber, double D_tot, double NRM);
void WriteStatisticsToFile(clock_t t, int orNumber, double D_tot, double NRM);

void TraceRandom(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				 Tracing &tracer);
void TraceFixed(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				Tracing &tracer);
void TraceSingle(Tracing &tracer, double beta, double gamma);

void ImportTracks(Particle *particle)
{
	const int bufSize = 1024;
	ifstream trackFile("tracks.dat", ios::in);
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	while (!trackFile.eof())
	{
		trackFile.getline(buff, bufSize);

		if (buff[0] == '@')
		{
			char *ptr = strtok(buff, "@ ");
			++groupCount;
			trackGroups[groupCount-1].groupID = strtol(ptr, &ptr, 10);
		}
		else
		{
			vector<int> arr;

			char *ptr, *trash;
			ptr = strtok(buff, " ");

			while (ptr != NULL)
			{
				int tmp = strtol(ptr, &trash, 10);
				arr.push_back(tmp);
				ptr = strtok(NULL, " ");
			}

			long long int trackID = 0;

			for (int t : arr)
			{
				trackID += (t + 1);
				trackID *= (particle->facetNum + 1);
			}

			trackGroups[groupCount-1].arr[trackGroups[groupCount-1].size++] = trackID;
			arr.clear();
		}
	}
}

void Calculate(const CLArguments &params)
{
	Tracing *tracer = nullptr;

	int orNumBeta = params.betaRange.end - params.betaRange.begin;
	int orNumGamma = params.gammaRange.end - params.gammaRange.begin;

	Particle *particle = nullptr;

	switch (params.particleType)
	{
	case ParticleType::Hexagonal:
		particle = new Hexagonal(params.radius, params.halfHeight, params.refractionIndex);
		tracer = new TracingConvex(particle, incidentDir, isOpticalPath,
								   polarizationBasis, params.interReflNum);
		betaNorm = M_PI/(2.0*orNumBeta);
		gammaNorm = M_PI/(3.0*orNumGamma);
		break;
	case ParticleType::TiltedHexagonal:
		particle = new TiltedHexagonal(params.radius, params.halfHeight, params.refractionIndex,
									   params.tiltAngle);
		tracer = new TracingConvex(particle, incidentDir, isOpticalPath,
								   polarizationBasis, params.interReflNum);
		betaNorm = M_PI/(2.0*orNumBeta);
		gammaNorm = (2.0*M_PI)/orNumGamma;
		break;
	case ParticleType::ConcaveHexagonal:
		particle = new ConcaveHexagonal(params.radius, params.halfHeight, params.refractionIndex,
										params.cavityDepth);
		tracer = new TracingConcave(particle, incidentDir, isOpticalPath,
									polarizationBasis, params.interReflNum);
		betaNorm = M_PI/(2.0*orNumBeta);
		gammaNorm = M_PI/(3.0*orNumGamma);
		break;
	default:
		assert(false && "ERROR! Incorrect type of particle.");
		break;
	}

	if (isPhisOptics)
	{
		ImportTracks(particle);

		double radius = (params.bsCone.radius*M_PI)/180.0;

		if (params.bsCone.theta)
		{
			dt = radius/(double)params.bsCone.theta;
		}

		if (params.bsCone.phi)
		{
			df = M_2PI/(double)(params.bsCone.phi+1);
		}
	}

	int EDF = 0;

	incomingEnergy = 0;

	// the arrays for exact backscattering and forwardscattering Mueller matrices
	back.Fill(0);
	forw.Fill(0);

	sizeBin = M_PI/params.thetaNumber; // the size of the bin (radians)

	mxd = Arr2D(1, params.thetaNumber+1, 4, 4);
	mxd.ClearArr();

	cout << endl;

	time_point<system_clock> startCalc = system_clock::now();

	if (params.isRandom)
	{
		TraceRandom(params.gammaRange, params.betaRange, *tracer);
	}
	else
	{
		TraceFixed(params.gammaRange, params.betaRange, *tracer);
	}

	auto total = system_clock::now() - startCalc;

	if (assertNum > 0)
	{
		cout << endl << "WARRNING! Asserts are occured (see log file) " << assertNum << endl;
	}

	cout << "\nTotal time of calculation = " << duration_cast<seconds>(total).count() << " seconds";

	// Integrating
	double D_tot = back[0][0] + forw[0][0];

	for (int j = 0; j <= params.thetaNumber; ++j)
	{
		D_tot += mxd(0, j, 0, 0);
	}

	// Normalizing coefficient
	long long orNum = orNumGamma * orNumBeta;
	double NRM;

	if (params.isRandom)
	{
		NRM = 1.0/(double)orNum;
	}
	else
	{
		NRM = M_PI/((double)orNum*2.0);
	}

	ExtractPeaks(EDF, NRM, params.thetaNumber);

	WriteResultsToFile(params.thetaNumber, NRM, params.outfile);
//	WriteStatisticsToFile(timer, orNum, D_tot, NRM);
	WriteStatisticsToConsole(orNum, D_tot, NRM);

	delete particle;
}

int GetArgValue(char* argv[], int argc, int i)
{
	if (argc <= i)
	{
		throw string("Not enouth arguments");
	}

	char *end;
	int val = strtol(argv[i], &end, 10);

	if (strlen(end) != 0)
	{
		throw string("Some argument is incorrect.");
	}

	return val;
}

double GetArgValueD(char* argv[], int argc, int i)
{
	if (argc <= i)
	{
		throw string("Not enouth arguments");
	}

	char *end;
	double val = strtod(argv[i], &end);

	if (strlen(end) != 0)
	{
		throw string("Some argument is incorrect.");
	}

	return val;
}

void SetParams(int argc, char* argv[], CLArguments &params)
{
	try
	{
		int paramsNum = 0;

		for (int i = 1; i < argc; ++i)
		{
			string arg(argv[i]);

			if (arg == "-p")
			{
				params.particleType = (ParticleType)GetArgValue(argv, argc, ++i);
				params.halfHeight = GetArgValueD(argv, argc, ++i);
				params.radius = GetArgValueD(argv, argc, ++i);

				if (params.particleType == ParticleType::ConcaveHexagonal)
				{
					params.cavityDepth = GetArgValueD(argv, argc, ++i);
				}
				else if (params.particleType == ParticleType::TiltedHexagonal)
				{
					params.tiltAngle = GetArgValueD(argv, argc, ++i);
				}

				++paramsNum;
			}
			else if (arg == "-ri")
			{
				params.refractionIndex = GetArgValueD(argv, argc, ++i);
				++paramsNum;
			}
			else if (arg == "-rn")
			{
				params.interReflNum = GetArgValue(argv, argc, ++i);
				++paramsNum;
			}
			else if (arg == "-b")
			{
				params.betaRange.begin = GetArgValue(argv, argc, ++i);
				params.betaRange.end = GetArgValue(argv, argc, ++i);
				++paramsNum;
			}
			else if (arg == "-g")
			{
				params.gammaRange.begin = GetArgValue(argv, argc, ++i);
				params.gammaRange.end = GetArgValue(argv, argc, ++i);
				++paramsNum;
			}
			else if (arg == "-t")
			{
				params.thetaNumber = GetArgValue(argv, argc, ++i);
				++paramsNum;
			}
			else if (arg == "-r")
			{
				params.isRandom = true;
			}
			else if (arg == "-o")
			{
				if (argc <= i)
				{
					throw "Not enouth arguments.";
				}

				params.outfile = argv[++i];
			}
			else if (arg == "-bsc")
			{
				params.bsCone.radius = GetArgValueD(argv, argc, ++i);
				params.bsCone.phi = GetArgValue(argv, argc, ++i);
				params.bsCone.theta = GetArgValue(argv, argc, ++i);
			}
			else if (arg == "-w")
			{
				params.wavelength = GetArgValueD(argv, argc, ++i);
			}
		}

		if (paramsNum < 6) // REF: выделить как константу
		{
			throw string("Too few arguments.");
		}
	}
	catch (const string &e)
	{
		cout << "Error! " << e << " Please check it and restart the program."
				  << endl << "Press any key to exit...";
		getchar();
		exit(1);
	}
}

int main(int argc, char* argv[])
{
//	logfile->open("log.txt", ios::out);

//	testConcaveHexagonRot();
//	testHexagonBuilding();
//	testHexagonRotate();
//	testTiltHexagonBuild();
//	testCompareParticles();

	if (argc > 1) // has command line arguments
	{
		SetParams(argc, argv, params);
	}
	else // DEB
	{
		cout << "Argument list is not found. Using default params."
				  << endl << endl;

//		params.particleType = ParticleType::TiltedHexagonal;
		params.particleType = ParticleType::Hexagonal;
		params.halfHeight = 100;
		params.radius = 40;
		params.cavityDepth = 20;
		params.tiltAngle = 45;
		params.refractionIndex = complex(1.31, 0.0);
		params.betaRange.begin = 0;
		params.betaRange.end = 200;
		params.gammaRange.begin = 0;
		params.gammaRange.end = 201;
		params.thetaNumber = 180;
		params.interReflNum = 3;
		params.isRandom = false;

#ifdef _OUTPUT_NRG_CONV
		cout << "WARNING: Energy conversation is calculating now."
				  << endl << endl;
		params.refractionIndex = complex(1000000000000001.31, 0.0);
		params.interReflNum = 10;
#endif
	}

	Calculate(params);
	getchar();
	return 0;
}

void TraceSingle(Tracing &tracer, double beta, double gamma)
{
	beta = (M_PI*beta)/180.0;
	gamma = (M_PI*gamma)/180.0;
	double betaDistrProbability = sin(beta);

	vector<Beam> outcomingBeams;
//	double square;

//	tracer.RotateParticle(beta, gamma);
	tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);

//	incomingEnergy += tracer.GetLightSurfaceArea();
	HandleBeams(outcomingBeams, betaDistrProbability, tracer);

	outcomingBeams.clear();
}

void PrintTime(long long &msLeft, CalcTimer &time)
{
	time.Left(msLeft);
	cout << "time left: " << time.ToString();
	cout << "\t\tends at " << ctime(&time.End(msLeft));
}

int GetMaxGroupID()
{
	int maxGroupID = 0;

	for (int i = 0; i < groupCount; ++i)
	{
		if (trackGroups[i].groupID > maxGroupID)
		{
			maxGroupID = trackGroups[i].groupID;
		}
	}

	return maxGroupID;
}

void WriteSumMatrix(ofstream &M_all_file, const Arr2D &M_)
{
	M_all_file << to_string(params.bsCone.radius) << ' '
			<< to_string(params.bsCone.theta) << ' '
			<< to_string(params.bsCone.phi+1);

	for (int t = 0; t <= params.bsCone.theta; ++t)
	{
		double tt = (double)(t*dt*180.0)/M_PI;

		for (int p = 0; p <= params.bsCone.phi; ++p)
		{
			double fi = -((double)p)*df;
			matrix m = M_(p ,t);
			M_all_file << endl << tt << " " << (-fi*180)/M_PI << " ";
			M_all_file << m;
		}
	}
}

void AddResultToSumMatrix(Arr2D &M_, int maxGroupID, double norm)
{
	for (int q = 0; q < maxGroupID; ++q)
	{
		for (int t = 0; t <= params.bsCone.theta; ++t)
		{
			for (int p = 0; p <= params.bsCone.phi; ++p)
			{
//				complex ee = J[q](p, t)[0][0];
				matrix Mk = Mueller(J[q](p, t));
//				M[q].insert(p, t, gammaNorm*norm*Mk);
				M_.insert(p, t, gammaNorm*norm*Mk);
			}
		}
	}
}

void CleanJ(int maxGroupID)
{
	J.clear();
	Arr2DC tmp(params.bsCone.phi+1, params.bsCone.theta+1, 2, 2);
	tmp.ClearArr();

	for(int q = 0; q < maxGroupID; q++)
	{
		J.push_back(tmp);
	}
}

void CleanM(vector<Arr2D> &M, int maxGroupID)
{
	M.clear();
	Arr2D M_void(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4);

	for (int j = 0; j < maxGroupID; ++j)
	{
		M.push_back(M_void);
	}
}

void WriteToSepatateFiles(vector<Arr2D> &M, double beta, int maxGroupID)
{
	for (unsigned int q = 0; q < maxGroupID; ++q)
	{
		string tr = "3673_3763";
		string orFName = ("b_" + to_string((beta*180.0)/M_PI) + "_" + tr + ".dat").c_str();
		ofstream or_file(orFName, ios::out);

		or_file << to_string(params.bsCone.radius) << ' '
				<< to_string(params.bsCone.theta ) << ' '
				<< to_string(params.bsCone.phi+1);

//				matrix sum(4, 4);

		for (int t = 0; t <= params.bsCone.theta; ++t)
		{
//					sum.Fill(0);
			double tt = (double)(t*dt*180.0)/M_PI;

			for (int p = 0; p <= params.bsCone.phi; ++p)
			{
				double fi = -((double)p)*df;
				matrix m = M[q](p ,t);

				or_file << endl << tt << " " << (-fi*180)/M_PI << " "; or_file << m;

//						matrix L(4,4);
//						L[0][0] = 1.0;
//						L[0][1] = 0.0;
//						L[0][2] = 0.0;
//						L[0][3] = 0.0;
//						L[1][0] = 0.0;
//						L[1][1] = cos(2.0*p);
//						L[1][2] = sin(2.0*p);
//						L[1][3] = 0.0;
//						L[2][0] = 0.0;
//						L[2][1] =-sin(2.0*p);
//						L[2][2] = cos(2.0*p);
//						L[2][3] = 0.0;
//						L[3][0] = 0.0;
//						L[3][1] = 0.0;
//						L[3][2] = 0.0;
//						L[3][3] = 1.0;

//						if (!t)
//						{
//							sum += L*m*L;
//						}
//						else
//						{
//							sum += m*L;
//						}
			}

//					sum /= (params.bsCone.phi+1.0);
//					M_all_file << endl << (beta*180.0)/M_PI << " " << tt << " ";
//					M_all_file << sum;
		}

		or_file.close();
//				M_all_file << endl;
	}
}

void TraceSinglePO(Tracing &tracer, double beta, double gamma)
{
	beta = (32*M_PI)/180;
	gamma = (30*M_PI)/180;

	vector<Arr2D> M;
	Arr2D M_(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4);

	int maxGroupID = GetMaxGroupID();
	double norm = 1.0/(M_PI/3.0);

	CleanM(M, maxGroupID);
	CleanJ(maxGroupID);

	ofstream M_all_file("M_all.dat", ios::out); // матрица Мюллера общая (физ. опт.)

	try
	{
		vector<Beam> outcomingBeams;
		tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);
		HandleBeams(outcomingBeams, 0, tracer);

		if (isPhisOptics)
		{
			AddResultToSumMatrix(M_, maxGroupID, norm);
		}
	}
	catch (const bool &)
	{
		logfile.flush();
		++assertNum;
	}

	if (isPhisOptics)
	{
		WriteSumMatrix(M_all_file, M_);
//		WriteToSepatateFiles(M, beta, maxGroupID);
	}

	M_all_file.close();
}

void TraceFixed(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				Tracing &tracer)
{
	CalcTimer time;

	double beta, gamma;
	double betaDistrProbability;

	vector<Beam> outcomingBeams;
	double square = 0;

	int orNumBeta = betaRange.end - betaRange.begin;
	int orNumGamma = gammaRange.end - gammaRange.begin;
	long long orNum = orNumGamma * orNumBeta;

	long long count = 0;

	double dBettaRad = 0.0;
	{
		if (orNumBeta)
		{
			dBettaRad = ((betaMax*M_PI)/180.0)/(double)orNumBeta;
		}
	}

	// PO params
	ofstream M_all_file("M_all.dat", ios::out); // матрица Мюллера общая (физ. опт.)
	int maxGroupID = GetMaxGroupID();
	++maxGroupID;

	vector<Arr2D> M;
	Arr2D M_(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4);
	//

	double norm = 1.0/(M_PI/3.0);

	time.Start();

#ifdef _DEBUG // DEB
{
	//	beta = (153 + 0.5)*betaNorm;
	//	gamma = (100 + 0.5)*gammaNorm;
	//	betaDistrProbability = sin(beta);
	//	tracer.RotateParticle(beta, gamma);
	//	tracer.SplitBeamByParticle(outcomingBeams);
	//	HandleBeams(outcomingBeams, betaDistrProbability, tracer);
}
#ifdef _OUTPUT_NRG_CONV
	double sss=3.0*sqrt(3.0)/2.0*40.0*40.0*sin(M_PI/2.0-beta)+2.0*40.0*200.0*cos(M_PI/2.0-beta)*cos(M_PI/6.0-gamma);
	energyFile<<153<<" "<<100<<" "<<beta*180.0/3.1415926<<" "<<gamma*180./3.1415926<<" "<<bcount<<" "<<sss<<" "<<SS<<" "<<(fabs(sss-SS)<0.1?0:sss-SS)<<endl;
	if (fabs(sss-SS) > 40)
		int fff = 0;
#endif
#endif

	double dbeta = (M_PI/2 - 0)/orNumBeta;

	for (int i = betaRange.begin; i <= betaRange.end; ++i)
	{
//		beta = (i + 0.5)*betaNorm;
		beta = 0 + (double)i*dbeta;
		betaDistrProbability = sin(beta);

		CleanM(M, maxGroupID);

		for (int j = gammaRange.begin; j <= gammaRange.end; ++j)
		{
//			gamma = (j + 0.5)*gammaNorm;
			gamma = (j - gammaRange.end/2)*gammaNorm + M_PI/6;

#ifdef _OUTPUT_NRG_CONV
			SS=0;
			bcount=0;
#endif
			CleanJ(maxGroupID);

			try
			{
//				tracer.RotateParticle(beta, gamma);
				tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);

				/// TODO: сделать отдельную ф-цию для расчёта фиксированных траекторий
	//			vector<vector<int>> tracks;
	//			vector<int> track = {0, 7, 0};
	//			tracks.push_back(track);
	//			tracer.SplitBeamByParticle(incidentDir, tracks, outcomingBeams);

				incomingEnergy += betaDistrProbability * square;
				HandleBeams(outcomingBeams, betaDistrProbability, tracer);

#ifdef _OUTPUT_NRG_CONV
				double sss=3.0*sqrt(3.0)/2.0*40.0*40.0*sin(M_PI/2.0-beta)+2.0*40.0*200.0*cos(M_PI/2.0-beta)*cos(M_PI/6.0-gamma);
				energyFile<<i<<" "<<j<<" "<<beta*180.0/3.1415926<<" "<<gamma*180./3.1415926<<" "<<bcount<<" "<<sss<<" "<<SS<<" "<<(fabs(sss-SS)<0.1?0:sss-SS)<<endl;
#ifdef _DEBUG // DEB
				if (fabs(sss-SS) > 40)
					int fff = 0;
#endif
#endif
				M_all_file << (beta*180)/M_PI << ' '
						   << (gamma*180)/M_PI << ' '
						   << ((gammaNorm*norm*Mueller(J[0](0, 0)))[0][0]) << endl;

				if (isPhisOptics)
				{
					AddResultToSumMatrix(M_, maxGroupID, norm);
				}
			}
			catch (const bool &)
			{
				logfile << ", i: " << i << ", j: " << j << endl;
				logfile.flush();
				++assertNum;
			}

			outcomingBeams.clear();
		}

		if (/*isPhisOptics*/false)
		{
			WriteSumMatrix(M_all_file, M_);
//			WriteToSepatateFiles(M, beta, maxGroupID);
		}

		Dellines(2);
		cout << ((100*count)/orNumBeta) << "% ";

		time.Stop();
		long long durMs = time.Duration();
		long long msLeft = (durMs/orNumGamma)*(orNum - count*orNumGamma);

		PrintTime(msLeft, time);

		time.Start();
		++count;
	}

	M_all_file.close();
}

void TraceRandom(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				 Tracing &tracer)
{
	srand(time(NULL));

	vector<Beam> outcomingBeams;
	double beta, gamma;
//	double square;
	int orNumBeta = betaRange.end - betaRange.begin;
	int count = 0;

	for (int i = betaRange.begin; i < betaRange.end; ++i)
	{
		for (int j = gammaRange.begin; j < gammaRange.end; ++j)
		{
			gamma = 2.0*M_PI * ((double)(rand()) / ((double)RAND_MAX + 1.0));
			beta = acos(((double)(rand()) / ((double)RAND_MAX))*2.0 - 1.0);

//			tracer.RotateParticle(beta, gamma);
			tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);

//			incomingEnergy += tracer.GetLightSurfaceArea();
			HandleBeams(outcomingBeams, 1, tracer);
			outcomingBeams.clear();
		}

		cout << (100*(++count))/orNumBeta << "%" << endl;
	}
}

void ExtractPeaks(int EDF, double NRM, int ThetaNumber)
{
	//Analytical averaging over alpha angle
	double b[3], f[3];
	b[0] = back[0][0];
	b[1] = (back[1][1] - back[2][2])/2.0;
	b[2] = back[3][3];

	f[0] = forw[0][0];
	f[1] = (forw[1][1] + forw[2][2])/2.0;
	f[2] = forw[3][3];

	// Extracting the forward and backward peak in a separate file if needed
	if (EDF)
	{
		ofstream bck("back.dat", ios::out);
		ofstream frw("forward.dat", ios::out);
		frw << "M11 M22/M11 M33/M11 M44/M11";
		bck << "M11 M22/M11 M33/M11 M44/M11";

		if (f[0] <= DBL_EPSILON)
		{
			frw << "\n0 0 0 0";
		}
		else
		{
			frw << "\n" << f[0]*NRM
				<< " " << f[1]/f[0]
				<< " " << f[1]/f[0]
				<< " " << f[2]/f[0];
		}

		if (b[0] <= DBL_EPSILON)
		{
			bck << "\n0 0 0 0";
		}
		else
		{
			bck << "\n" << b[0]*NRM
				<< " " << b[1]/b[0]
				<< " " << -b[1]/b[0]
				<< " " << b[2]/b[0];
		}

		bck.close();
		frw.close();
	}
	else
	{
		mxd(0,ThetaNumber,0,0) += f[0];
		mxd(0,0,0,0) += b[0];
		mxd(0,ThetaNumber,1,1) += f[1];
		mxd(0,0,1,1) += b[1];
		mxd(0,ThetaNumber,2,2) += f[1];
		mxd(0,0,2,2) -= b[1];
		mxd(0,ThetaNumber,3,3) += f[2];
		mxd(0,0,3,3) += b[2];
	}
}

bool IsMatchTrack(const vector<int> &track, const vector<int> &compared)
{
	if (track.size() != compared.size())
	{
		return false;
	}

	for (unsigned int i = 0; i < compared.size(); ++i)
	{
		if (track.at(i) != compared.at(i))
		{
			return false;
		}
	}

	return true;
}

int GetGroupID(long long int track)
{
	for (int i = 0; i < groupCount; ++i)
	{
		for (int j = 0; j < trackGroups[i].size; ++j)
		{
			if (trackGroups[i].arr[j] == track)
			{
				return trackGroups[i].groupID;
			}
		}
	}

	return -1;
}

void HandleBeams(vector<Beam> &outBeams, double betaDistrProb, const Tracing &tracer)
{
#ifdef _DEBUG // DEB
//	double eee = 0;
#endif

	if (isPhisOptics)
	{
		for (unsigned int i = 0; i < outBeams.size(); ++i)
		{
			Beam &beam = outBeams.at(i);

#ifdef _DEBUG // DEB
//			eee += beam.polygon.Area();
#endif
			double ctetta = DotProduct(beam.direction, -incidentDir);

			if (ctetta < 0.17364817766693034885171662676931)
			{	// отбрасываем пучки, которые далеко от конуса направления назад
//				continue;// DEB
			}

			int groupID = GetGroupID(beam.id);

			if (groupID < 0)
			{
				continue;
			}

			beam.RotateSpherical(-incidentDir, polarizationBasis);

			Point3f center = beam.polygon.Center();
			double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

			Point3f T = CrossProduct(beam.e, beam.direction);
			T = T/Length(T); // базис выходящего пучка

			ofstream file("dd.dat", ios::out); // DEB
			file << beam.id << endl;

			for (int i = 0; i <= params.bsCone.phi; ++i)
			{
				for (int j = 0; j <= params.bsCone.theta; ++j)
				{
					double f = (double)i*df;
					double t = (double)j*dt;

					double sinT = sin(t),
							sinF = sin(f),
							cosF = cos(f);

					Point3d vr(sinT*cosF, sinT*sinF, cos(t));

					matrixC Jn_rot(2, 2);
					{
						Point3f normal = beam.polygon.Normal();

						Point3d vf, vt;
						vf = (!j) ? -polarizationBasis : Point3d(-sinF ,cosF ,0);
						vt = CrossProductD(vf, vr);
						vt = vt/LengthD(vt);

						Point3f cpNT = CrossProduct(normal, T);
						Point3f cpNE = CrossProduct(normal, beam.e);

						Jn_rot[0][0] = -DotProductD(Point3d(cpNT.cx, cpNT.cy, cpNT.cz), vf); // OPT: похоже на SetJMatrix
						Jn_rot[0][1] = -DotProductD(Point3d(cpNE.cx, cpNE.cy, cpNE.cz), vf);
						Jn_rot[1][0] =  DotProductD(Point3d(cpNT.cx, cpNT.cy, cpNT.cz), vt);
						Jn_rot[1][1] =  DotProductD(Point3d(cpNE.cx, cpNE.cy, cpNE.cz), vt);
					}

					complex fn(0, 0);

					// DEB
//					if (i == 90 && j == 15 /*&& beam.id == 32679*/)
//						int fff = 0;

					fn = beam.DiffractionIncline(vr, params.wavelength);

					double dp = DotProductD(vr, Point3d(center.cx, center.cy, center.cz));
					complex tmp = exp_im(M_2PI*(lng_proj0-dp)/params.wavelength);
					matrixC fn_jn = beam.J * tmp;

					matrixC c = fn*Jn_rot*fn_jn;
//					complex d1 = fn_jn[0][0];
//					complex d2 = fn_jn[0][1];
//					complex d3 = fn_jn[1][0];
//					complex d4 = fn_jn[1][1];
					J[groupID].insert(i, j, fn*Jn_rot*fn_jn);

					if (i == 0 && j == 0)
					{
						complex ff = c[0][0];
						file << i << ' ' << j << ' ' << real(ff) << ' ' << imag(ff)<< endl;
					}
				}
			}
		}
	}
	else
	{
		for (unsigned int i = 0; i < outBeams.size(); ++i)
		{
			Beam &beam = outBeams.at(i);

			beam.RotateSpherical(incidentDir, polarizationBasis);

			double cross = tracer.BeamCrossSection(beam);
			double Area = betaDistrProb * cross;

#ifdef _DEBUG // DEB
//			eee += Area;
#endif
			matrix bf = Mueller(beam.J);

			const float &x = beam.direction.cx;
			const float &y = beam.direction.cy;
			const float &z = beam.direction.cz;

#ifdef _OUTPUT_NRG_CONV
			if (bf[0][0]>0.000001)
			{
				bcount++;
				SS+=cross;
			}
#endif
			// Collect the beam in array
			if (z >= 1-DBL_EPSILON)
			{
				back += Area*bf;
			}
			else if (z <= DBL_EPSILON-1)
			{
				forw += Area*bf;
			}
			else
			{
				// Rotate the Mueller matrix of the beam to appropriate coordinate system
				const unsigned int ZenAng = round(acos(z)/sizeBin);
				double tmp = y*y;

				if (tmp > DBL_EPSILON)
				{
					tmp = acos(x/sqrt(x*x+tmp));

					if (y < 0)
					{
						tmp = M_2PI-tmp;
					}

					tmp *= -2.0;
					RightRotateMueller(bf, cos(tmp), sin(tmp));
				}

#ifdef _CALC_AREA_CONTIBUTION_ONLY
				bf = matrix(4,4);
				bf.Identity();
#endif
				//			matrix cc = Area*bf;
				mxd.insert(0, ZenAng, /*cc*/Area*bf);
			}

			LOG_ASSERT(Area >= 0);
		}
	}

#ifdef _DEBUG // DEB
//		int fff = 0;
//		complex ff = (J[0](0, 0))[0][0];
//		cout << endl << real(ff) << ' ' << imag(ff) << endl << endl << endl << endl;
#endif
}

string GetFileName(const string &filename)
{
	string fname = string("M_") + filename;
	string name = fname;

	for (int i = 1; ifstream(name += ".dat") != NULL; ++i)
	{
		name = fname + "(" + to_string(i) + ")";
	}

	return name;
}

void WriteResultsToFile(int ThetaNumber, double NRM, const string &filename)
{
	string name = GetFileName(filename);

	ofstream M(name, ios::out);

	M <<  "tetta M11 M12/M11 M21/M11 M22/M11 M33/M11 M34/M11 M43/M11 M44/M11";

	for (int j = ThetaNumber; j >= 0; j--)
	{
		double sn;

		//Special case in first and last step
		M << '\n' << 180.0/ThetaNumber*(ThetaNumber-j) + (j==0 ?-0.25*180.0/ThetaNumber:0)+(j==(int)ThetaNumber ?0.25*180.0/ThetaNumber:0);
		sn = (j==0 || j==(int)ThetaNumber) ? 1-cos(sizeBin/2.0) : (cos((j-0.5)*sizeBin)-cos((j+0.5)*sizeBin));

		matrix bf = mxd(0,j);

		if(bf[0][0] <= DBL_EPSILON)
		{
			M << " 0 0 0 0 0 0 0 0";
		}
		else
		{
			M << ' ' << bf[0][0]*NRM/(2.0*M_PI*sn)
					<< ' ' << bf[0][1]/bf[0][0]
					<< ' ' << bf[1][0]/bf[0][0]
					<< ' ' << bf[1][1]/bf[0][0]
					<< ' ' << bf[2][2]/bf[0][0]
					<< ' ' << bf[2][3]/bf[0][0]
					<< ' ' << bf[3][2]/bf[0][0]
					<< ' ' << bf[3][3]/bf[0][0];
		}
	}
	M.close();
}

void WriteStatisticsToFile(clock_t t, int orNumber, double D_tot, double NRM)
{
	ofstream out("out.dat", ios::out);
	// Information for log-file
	out << "\nTotal time of calculation = " << t/CLOCKS_PER_SEC << " seconds";
	out << "\nTotal number of body orientation = " << orNumber;
	out << "\nTotal scattering energy = " << D_tot;
	out << "\nTotal incoming energy = " << incomingEnergy;
	out << "\nAveraged cross section = " << incomingEnergy*NRM;
	out.close();
}

void WriteStatisticsToConsole(int orNumber, double D_tot, double NRM)
{
	using namespace std;
	cout << "\nTotal number of body orientation = " << orNumber;
	cout << "\nTotal scattering energy = " << D_tot/**NRM*/;
	cout << "\nTotal incoming energy = " << incomingEnergy;
	cout << "\nAveraged cross section = " << incomingEnergy*NRM;
	cout << "\nAll done. Please, press ENTER.";
}
