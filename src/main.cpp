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

void ImportTracks2(int facetNum)
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

//void ImportTracks(int facetNum)
//{
//	const int bufSize = 1024;
//	ifstream trackFile("tracks.dat", ios::in);
//	char *buff = (char*)malloc(sizeof(char) * bufSize);

//	while (!trackFile.eof())
//	{
//		trackFile.getline(buff, bufSize);

//		if (buff[0] == '@')
//		{
//			char *ptr = strtok(buff, "@ ");
//			++trackGroups.count;
//			trackGroups.groups[trackGroups.count-1].groupID = strtol(ptr, &ptr, 10);
//		}
//		else
//		{
//			vector<int> track;

//			char *ptr, *trash;
//			ptr = strtok(buff, " ");

//			while (ptr != NULL)
//			{
//				int tmp = strtol(ptr, &trash, 10);
//				track.push_back(tmp);
//				ptr = strtok(NULL, " ");
//			}

//			trackGroups.groups[trackGroups.count-1].tracks.push_back(track);

//			long long int trackID = 0;

//			for (int t : track)
//			{
//				trackID += (t + 1);
//				trackID *= (facetNum + 1);
//			}

//			trackGroups.groups[trackGroups.count-1].arr[trackGroups.groups[trackGroups.count-1].size++] = trackID;
//			track.clear();
//		}
//	}
//}

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
//	logfile->open("log.txt", ios::out);

//	testConcaveHexagonRot();
//	testHexagonBuilding();
//	testHexagonRotate();
//	testTiltHexagonBuild();
//	testCompareParticles();

//	testHexagonalAggregateBuild();
//	testHexagonalAggregateRot(0, 0);
//	testHexagonalAggregateRot(45, 90);
//	testHexagonalAggregateRot(45, -90);
//	testHexagonalAggregateRot(30/*RadToDeg(0.001963495408493621)*/,
//							  RadToDeg(0.00392208820672883));

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
		ImportTracks2(particle->facetNum);
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
					tracer.TraceBackScatterPointPO(betaR, gammaR, trackGroups, 0.532);
//					tracer.TraceIntervalGO(betaR, gammaR, cellNum, trackGroups);
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

void WriteSumMatrix(ofstream &M_all_file, const Arr2D &M_)
{
//	M_all_file << to_string(params.bsCone.radius) << ' '
//			<< to_string(params.bsCone.theta) << ' '
//			<< to_string(params.bsCone.phi+1);

	for (int t = 0; t <= 0 /*params.bsCone.theta*/; ++t)
	{
		double tt = (double)(t*dt*180.0)/M_PI;

		for (int p = 0; p <= 0/*params.bsCone.phi*/; ++p)
		{
			double fi = -((double)p)*df;
			matrix m = M_(p ,t);
			M_all_file << endl << tt << " " << (-fi*180)/M_PI << " ";
			M_all_file << m;
		}
	}
}

void AddResultToSumMatrix(Arr2D &M, int maxGroupID, double norm)
{
	for (int q = 0; q < maxGroupID; ++q)
	{
		for (int t = 0; t <= 0/*params.bsCone.theta*/; ++t)
		{
			for (int p = 0; p <= 0/*params.bsCone.phi*/; ++p)
			{
//				complex ee = J[q](p, t)[0][0];
				matrix Mk = Mueller(J[q](p, t));
//				M[q].insert(p, t, gammaNorm*norm*Mk);
				M.insert(p, t, gammaNorm*norm*Mk);
			}
		}
	}
}

void CleanM(vector<Arr2D> &M, int maxGroupID)
{
	M.clear();
	Arr2D M_void/*(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4)*/;

	for (int j = 0; j < maxGroupID; ++j)
	{
		M.push_back(M_void);
	}
}

void WriteToSepatateFiles(vector<Arr2D> &M, double beta, int maxGroupID)
{
	for (int q = 0; q < maxGroupID; ++q)
	{
		string tr = to_string(q);
		string orFName = (tr + "_" + "b_" + to_string((beta*180.0)/M_PI) + ".dat").c_str();
		ofstream or_file(orFName, ios::app);

//		or_file << to_string(params.bsCone.radius) << ' '
//				<< to_string(params.bsCone.theta ) << ' '
//				<< to_string(params.bsCone.phi+1);

//				matrix sum(4, 4);

		for (int t = 0; t <= 0/*params.bsCone.theta*/; ++t)
		{
//					sum.Fill(0);
			double tt = (double)(t*dt*180.0)/M_PI;

			for (int p = 0; p <= 0/*params.bsCone.phi*/; ++p)
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
	Arr2D M_/*(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4)*/;

	int maxGroupID /*= GetMaxGroupID()*/;
	double norm = 1.0/(M_PI/3.0);

	CleanM(M, maxGroupID);
//	CleanJ(maxGroupID);

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

	// PO params
	ofstream M_all_file("M_all.dat", ios::out); // матрица Мюллера общая (физ. опт.)
	int maxGroupID = /*GetMaxGroupID()*/0;
	++maxGroupID;

	vector<Arr2D> M;
//	Arr2D M_(params.bsCone.phi+1, params.bsCone.theta+1, 4, 4);
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

	for (int i = betaRange.begin; i <= betaRange.end/*betaRange.begin*/; ++i)
	{
//		beta = (i + 0.5)*betaNorm;
		beta = 0 + (double)i*dbeta;
		betaDistrProbability = sin(beta);

		CleanM(M, maxGroupID);

		double ddd = 0;

		for (int j = gammaRange.begin; j <= gammaRange.end/*gammaRange.begin*/; ++j)
		{
//			gamma = (j + 0.5)*gammaNorm;
			gamma = (j - gammaRange.end/2)*gammaNorm;

#ifdef _OUTPUT_NRG_CONV
			SS=0;
			bcount=0;
#endif
//			CleanJ(maxGroupID);

			try
			{
				tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);
//gamma = beta = 0;
//beta = (0.25*M_PI)/180;
//gamma = (50*M_PI)/180;

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
				for (int q = 0; q < maxGroupID; ++q)
				{
					ddd += ((gammaNorm*norm*Mueller(J[q](0, 0)))[0][0]);
				}

//				M_all_file << (beta*180)/M_PI << ' '
//						   << (gamma*180)/M_PI << ' '
//						   << ((gammaNorm*norm*Mueller(J[0](0, 0)))[0][0]) << endl;

				if (isPhisOptics)
				{
//					for (int k = 0; k < maxGroupID; ++k)
//					{
//						for (int t = 0; t <= params.bsCone.theta; ++t)
//						{
//							for (int p = 0; p <= params.bsCone.phi; ++p)
//							{
////								complex ee = J[k](p, t)[0][0];
//								matrix Mk = Mueller(J[k](p, t));
//								M[k].insert(p, t, gammaNorm*norm*Mk);
//							}
//						}
//					}
//					AddResultToSumMatrix(M_, maxGroupID, norm);
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

		M_all_file << (beta*180)/M_PI << ' '
				   << (gamma*180)/M_PI << ' '
				   << ddd << endl;

		if (/*isPhisOptics*/true)
		{
//			WriteSumMatrix(M_all_file, M_);
			WriteToSepatateFiles(M, beta, maxGroupID);
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

//void TraceRandom(const AngleInterval &gammaRange, const AngleInterval &betaRange,
//				 Tracing &tracer)
//{
//	srand(time(NULL));

//	vector<Beam> outcomingBeams;
//	double beta, gamma;
////	double square;
//	int orNumBeta = betaRange.end - betaRange.begin;
//	int count = 0;

//	for (int i = betaRange.begin; i < betaRange.end; ++i)
//	{
//		for (int j = gammaRange.begin; j < gammaRange.end; ++j)
//		{
//			gamma = 2.0*M_PI * ((double)(rand()) / ((double)RAND_MAX + 1.0));
//			beta = acos(((double)(rand()) / ((double)RAND_MAX))*2.0 - 1.0);

////			tracer.RotateParticle(beta, gamma);
//			tracer.SplitBeamByParticle(beta, gamma, outcomingBeams);

////			incomingEnergy += tracer.GetLightSurfaceArea();
//			HandleBeams(outcomingBeams, 1, tracer);
//			outcomingBeams.clear();
//		}

//		cout << (100*(++count))/orNumBeta << "%" << endl;
//	}
//}

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

			int groupID/* = GetGroupID(beam.id)*/;

			if (groupID < 0)
			{
				continue;
			}

			beam.RotateSpherical(-incidentDir, polarizationBasis);

			Point3f center = beam.Center();
			double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

			Point3f T = CrossProduct(beam.e, beam.direction);
			T = T/Length(T); // базис выходящего пучка

			ofstream file("dd.dat", ios::out); // DEB
			file << beam.id << endl;

			for (int i = 0; i <= 0/*params.bsCone.phi*/; ++i)
			{
				for (int j = 0; j <= 0/*params.bsCone.theta*/; ++j)
				{
					double f = (double)i*df;
					double t = (double)j*dt;

					double sinT = sin(t),
							sinF = sin(f),
							cosF = cos(f);

					Point3d vr(sinT*cosF, sinT*sinF, cos(t));

					matrixC Jn_rot(2, 2);
					{
						Point3f normal = beam.Normal();

						Point3d vf, vt;
						vf = (!j) ? -polarizationBasis : Point3d(-sinF ,cosF ,0);
						vt = CrossProductD(vf, vr);
						vt = vt/LengthD(vt);

						Point3f cpNT = CrossProduct(normal, T);
						Point3f cpNE = CrossProduct(normal, beam.e);

						Jn_rot[0][0] = -DotProductD(Point3d(cpNT.cx, cpNT.cy, cpNT.cz), vf);
						Jn_rot[0][1] = -DotProductD(Point3d(cpNE.cx, cpNE.cy, cpNE.cz), vf);
						Jn_rot[1][0] =  DotProductD(Point3d(cpNT.cx, cpNT.cy, cpNT.cz), vt);
						Jn_rot[1][1] =  DotProductD(Point3d(cpNE.cx, cpNE.cy, cpNE.cz), vt);
					}

					complex fn(0, 0);

					// DEB
//					if (i == 90 && j == 15 /*&& beam.id == 32679*/)
//						int fff = 0;

					fn = beam.DiffractionIncline(vr, 0/*params.wavelength*/);

					double dp = DotProductD(vr, Point3d(center.cx, center.cy, center.cz));
					complex tmp = exp_im(M_2PI*(lng_proj0-dp)/*/params.wavelength*/);
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

void WriteStatisticsToConsole(int orNumber, double D_tot, double NRM)
{
	using namespace std;
	cout << "\nTotal number of body orientation = " << orNumber;
	cout << "\nTotal scattering energy = " << D_tot/**NRM*/;
	cout << "\nTotal incoming energy = " << incomingEnergy;
	cout << "\nAveraged cross section = " << incomingEnergy*NRM;
	cout << "\nAll done. Please, press ENTER.";
}
