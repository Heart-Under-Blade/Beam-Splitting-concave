#include <iostream>
#include <time.h>

#include "test.h"

#include <assert.h>
#include <float.h>
#include "global.h"

#include "Mueller.hpp"

#include "Hexagonal.h"
#include "ConcaveHexagonal.h"

#include "TracingConcave.h"
#include "TracingConvex.h"

#include "Beam.h"
#include "PhysMtr.hpp"

struct OrientationRange
{
	int begin;
	int end;
};

enum class ParticleType : int
{
	Hexagonal = 1,
	ConcaveHexagonal = 10
};

struct CLArguments
{
	ParticleType particleType;
	double halfHeight;
	double radius;
	double cavityDepth;
	complex refractionIndex;
	OrientationRange betaRange;
	OrientationRange gammaRange;
	int thetaNumber;
	int interReflNum;
	bool isRandom = false;
	std::string outfile;
};

matrix back(4,4),	///< Mueller matrix in backward direction
		forw(4,4);	///< Mueller matrix in forward direction

//std::ofstream WW("WW.dat", std::ios::out); //DEB

Arr2D mxd(0, 0, 0, 0);
double SizeBin;
double gammaNorm, betaNorm;
double incomingEnergy;
clock_t totalTime;
Point3f incidentDir(0, 0, 1); /// REF: исправить на нормальное направление
Point3f polarizationBasis(0, 1, 0);

// DEB
long long ddddd = 0;

void HandleBeams(std::vector<Beam> &outBeams, double betaDistrProb, const Tracing &tracer);
void ExtractPeaks(int EDF, double NRM, int ThetaNumber);

void WriteResultsToFile(int ThetaNumber, double NRM);
void WriteStatisticsToConsole(int orNumber, double D_tot, double NRM);
void WriteStatisticsToFile(clock_t t, int orNumber, double D_tot, double NRM);

void TraceRandom(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				 Tracing &tracer);
void TraceFixed(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				Tracing &tracer);
void TraceSingle(Tracing &tracer, double beta, double gamma);

void Calculate(const CLArguments &params)
{
	bool isOpticalPath = false;

	Tracing *tracer = nullptr;
	Particle *particle = nullptr;

	int orNumBeta = params.betaRange.end - params.betaRange.begin;
	int orNumGamma = params.gammaRange.end - params.gammaRange.begin;

	switch (params.particleType) {
	case ParticleType::Hexagonal:
		particle = new Hexagonal(params.radius, params.halfHeight, params.refractionIndex);
		tracer = new TracingConvex(particle, incidentDir, isOpticalPath,
								   polarizationBasis, params.interReflNum);
		betaNorm = M_PI/(2.0*orNumBeta);
		break;
//	case 1:
		/// TODO реализовать остальные частицы
//		break;
	case ParticleType::ConcaveHexagonal:
//		particle = new Hexagonal(radius, halfHeight, refractionIndex); // DEB
		particle = new ConcaveHexagonal(params.radius, params.halfHeight,
										params.refractionIndex, params.cavityDepth);
		tracer = new TracingConcave(particle, incidentDir, isOpticalPath,
									polarizationBasis, params.interReflNum);
		betaNorm = M_PI/(2.0*orNumBeta); /// TODO: какое д/б betaNorm?
		break;
	default:
		break;
	}

	int EDF = 0;
	gammaNorm = M_PI/(3.0*orNumGamma);

	incomingEnergy = 0;

	// the arrays for exact backscattering and forwardscattering Mueller matrices
	back.Fill(0);
	forw.Fill(0);

	SizeBin = M_PI/params.thetaNumber; // the size of the bin (radians)

	mxd = Arr2D(1, params.thetaNumber+1, 4, 4);
	mxd.ClearArr();

	clock_t timer = clock();

	if (params.isRandom)
	{
		TraceRandom(params.gammaRange, params.betaRange, *tracer);
	}
	else
	{
		TraceFixed(params.gammaRange, params.betaRange, *tracer);
	}

	timer = clock() - timer;
	totalTime = timer/CLOCKS_PER_SEC;

	std::cout << "\nTotal time of calculation = " << totalTime << " seconds";

	// Integrating
	double D_tot = back[0][0] + forw[0][0];

	for (int j = 0; j <= params.thetaNumber; ++j)
	{
		D_tot += mxd(0, j, 0, 0);
	}

	// Normalizing coefficient
	double orNum = orNumGamma * orNumBeta;
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

	WriteResultsToFile(params.thetaNumber, NRM);
	WriteStatisticsToFile(timer, orNum, D_tot, NRM);
	WriteStatisticsToConsole(orNum, D_tot, NRM);

	delete particle;
}

int GetNextArgValue(char* argv[], int argc, int i)
{
	if (argc <= i)
	{
		throw std::string("Not enouth arguments");
	}

	char *end;
	int val = strtol(argv[i], &end, 10);

	if (strlen(end) != 0)
	{
		throw std::string("Some argument is incorrect.");
	}

	return val;
}

double GetArgValueD(char* argv[], int argc, int i)
{
	if (argc <= i)
	{
		throw std::string("Not enouth arguments");
	}

	char *end;
	double val = strtod(argv[i], &end);

	if (strlen(end) != 0)
	{
		throw std::string("Some argument is incorrect.");
	}

	return val;
}

int main(int argc, char* argv[])
{
//	testConcaveHexagonRot();
//	testHexagonBuilding();
//	testHexagonRotate();
// TODO: добавить файл лога и записывать туда ассерты
	CLArguments params;

	if (argc > 1) // has command line arguments
	{
		try
		{
			int paramsNum = 0;

			for (int i = 1; i < argc; ++i)
			{
				std::string arg(argv[i]);

				if (arg == "-p")
				{
					params.particleType = (ParticleType)GetNextArgValue(argv, argc, ++i);
					params.halfHeight = GetArgValueD(argv, argc, ++i);
					params.radius = GetArgValueD(argv, argc, ++i);

					if (params.particleType == ParticleType::ConcaveHexagonal)
					{
						params.cavityDepth = GetArgValueD(argv, argc, ++i);
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
					params.interReflNum = GetNextArgValue(argv, argc, ++i);
					++paramsNum;
				}
				else if (arg == "-b")
				{
					params.betaRange.begin = GetNextArgValue(argv, argc, ++i);
					params.betaRange.end = GetNextArgValue(argv, argc, ++i);
					++paramsNum;
				}
				else if (arg == "-g")
				{
					params.gammaRange.begin = GetNextArgValue(argv, argc, ++i);
					params.gammaRange.end = GetNextArgValue(argv, argc, ++i);
					++paramsNum;
				}
				else if (arg == "-t")
				{
					params.thetaNumber = GetNextArgValue(argv, argc, ++i);
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
			}

			if (paramsNum == 6) // REF: выделить как константу
			{
				throw std::string("Too few arguments.");
			}
		}
		catch (const std::string &e)
		{
			std::cout << "Error! " << e << " Please check it and restart the program."
					  << std::endl << "Press any key to exit...";
			getchar();
			exit(1);
		}

//		particleType = atoi(argv[i++]);
//		halfHeight = atof(argv[i++]);
//		radius = atof(argv[i++]);

//		if (particleType == 10)
//		{
//			cavityDept = atof(argv[i++]);
//		}

//		double re = atof(argv[i++]);
//		refractionIndex = complex(re, 0.0);

//		orNumber.gamma = atoi(argv[i++]);
//		orNumber.beta = atoi(argv[i++]);
//		thetaNumber = atoi(argv[i++]);
//		interReflNum = atoi(argv[i++]);
//		isRandom = atoi(argv[i++]);
	}
	else // DEB
	{
		std::cout << "Argument list is not found. Using default params."
				  << std::endl << std::endl;

		params.particleType = ParticleType::ConcaveHexagonal;
		params.halfHeight = 100;
		params.radius = 40;
		params.cavityDepth = 10;
		params.refractionIndex = complex(1.31, 0.0);
		params.betaRange.begin = 0;
		params.betaRange.end = 100;
		params.gammaRange.begin = 0;
		params.gammaRange.end = 101;
		params.thetaNumber = 180;
		params.interReflNum = 4;
		params.isRandom = false;
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

	std::vector<Beam> outcomingBeams;
	double square;

	tracer.RotateParticle(beta, gamma);
	tracer.SplitBeamByParticle(outcomingBeams, square);

	incomingEnergy += betaDistrProbability * square;
	HandleBeams(outcomingBeams, betaDistrProbability, tracer);

	outcomingBeams.clear();
}

void TraceFixed(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				Tracing &tracer)
{
	double beta, gamma;
	double betaDistrProbability;

	std::vector<Beam> outcomingBeams;
	double square = 0;

	// DEB
//	beta = (45 + 0.5)*betaNorm;
//	gamma = (54 + 0.5)*gammaNorm;
//	tracer.RotateParticle(beta, gamma);
//	tracer.SplitBeamByParticle(outcomingBeams, square);
//	HandleBeams(outcomingBeams, sin(beta), tracer);

	int orNumBeta = betaRange.end - betaRange.begin;
	int count = 0;

	for (int i = betaRange.begin; i < betaRange.end; ++i)
	{
		beta = (i + 0.5)*betaNorm;
		betaDistrProbability = sin(beta);

		for (int j = gammaRange.begin; j < gammaRange.end; ++j)
		{
#ifdef _DEBUG
//			std::cout << j << std::endl;
#endif
			gamma = (j + 0.5)*gammaNorm;

			tracer.RotateParticle(beta, gamma);
			tracer.SplitBeamByParticle(outcomingBeams, square);

			/// TODO: сделать отдельную ф-цию для расчёта фиксированных траекторий
//			std::vector<std::vector<int>> tracks;
//			std::vector<int> track = {0, 7, 0};
//			tracks.push_back(track);
//			tracer.SplitBeamByParticle(incidentDir, tracks, outcomingBeams);

			incomingEnergy += betaDistrProbability * square;
			HandleBeams(outcomingBeams, betaDistrProbability, tracer);
			outcomingBeams.clear();
		}

		std::cout << (100*(++count))/orNumBeta << "%" << std::endl;
	}
}

void TraceRandom(const OrientationRange &gammaRange, const OrientationRange &betaRange,
				 Tracing &tracer)
{
	srand(time(NULL));

	std::vector<Beam> outcomingBeams;
	double beta, gamma;
	double square;
	int orNumBeta = betaRange.end - betaRange.begin;
	int count = 0;

	for (int i = betaRange.begin; i < betaRange.end; ++i)
	{
		for (int j = gammaRange.begin; j < gammaRange.end; ++j)
		{
			gamma = 2.0*M_PI * ((double)(rand()) / ((double)RAND_MAX + 1.0));
			beta = acos(((double)(rand()) / ((double)RAND_MAX))*2.0 - 1.0);

			tracer.RotateParticle(beta, gamma);
			tracer.SplitBeamByParticle(outcomingBeams, square);

			incomingEnergy += square;
			HandleBeams(outcomingBeams, 1, tracer);
			outcomingBeams.clear();
		}

		std::cout << (100*(++count))/orNumBeta << "%" << std::endl;
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
		std::ofstream bck("back.dat", std::ios::out);
		std::ofstream frw("forward.dat", std::ios::out);
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

bool IsMatchTrack(const std::vector<int> &track, const std::vector<int> &compared)
{
	if (track.size() != compared.size())
	{
		return false;
	}

	for (int i = 0; i < compared.size(); ++i)
	{
		if (track.at(i) != compared.at(i))
		{
			return false;
		}
	}

	return true;
}

void HandleBeams(std::vector<Beam> &outBeams, double betaDistrProb, const Tracing &tracer)
{
	ddddd += outBeams.size();//DEB

	double ee = 0;//DEB
	for (unsigned int i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams.at(i);

		// DEB
//		if (!(IsMatchTrack(beam.track, {0})
//			  || IsMatchTrack(beam.track, {1})
//			|| IsMatchTrack(beam.track, {2})
//			  || IsMatchTrack(beam.track, {3})
//			  || IsMatchTrack(beam.track, {4})
//			|| IsMatchTrack(beam.track, {5})
//			  || IsMatchTrack(beam.track, {7})
//			  || IsMatchTrack(beam.track, {8})
//			|| IsMatchTrack(beam.track, {9})
//			  || IsMatchTrack(beam.track, {10})
//			|| IsMatchTrack(beam.track, {11})
//				/*|| IsMatchTrack(beam.track, {9,11,17,6,7})*/))
//		{
//			continue;
//		}

		beam.RotateSpherical(incidentDir, polarizationBasis);

		double cross = tracer.BeamCrossSection(beam);

		double Area = betaDistrProb * cross;
//		double Area = 1;// DEB
ee += Area;//DEB

		matrix bf = Mueller(beam.JMatrix);

		const float &x = beam.direction.cx;
		const float &y = beam.direction.cy;
		const float &z = beam.direction.cz;

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
			const unsigned int ZenAng = round(acos(z)/SizeBin);
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

			// DEB
//			bf = matrix(4,4);
//			bf.Identity();
//			if (ZenAng == 3 && (int)Area == 107)
			mxd.insert(0, ZenAng, Area*bf);

//			if (ZenAng == 3)
//				std::cout << (int)Area << std::endl;
//			count += (int)Area;
//			if (ZenAng == 3 /*&& Area > 941.782 && Area < 941.784*/)
//				WW << count << std::endl;//int fff = 0;

//			if (count < 0)
//				int fff = 0;
//			if (isinf(mxd(0, ZenAng, 0, 0)))
//				int fff = 0; // DEB
		}

		assert(Area >= 0);
	}
	ee = 0;// DEB
}

void WriteResultsToFile(int ThetaNumber, double NRM)
{
	std::ofstream M("M.dat", std::ios::out);

	M <<  "tetta M11 M12/M11 M21/M11 M22/M11 M33/M11 M34/M11 M43/M11 M44/M11";

	for (int j = ThetaNumber; j >= 0; j--)
	{
		double sn;

		//Special case in first and last step
		M << '\n' << 180.0/ThetaNumber*(ThetaNumber-j) + (j==0 ?-0.25*180.0/ThetaNumber:0)+(j==(int)ThetaNumber ?0.25*180.0/ThetaNumber:0);
		sn = (j==0 || j==(int)ThetaNumber) ? 1-cos(SizeBin/2.0) : (cos((j-0.5)*SizeBin)-cos((j+0.5)*SizeBin));

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
	std::ofstream out("out.dat", std::ios::out);
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
	cout << endl << ddddd << endl;
	cout << "\nTotal number of body orientation = " << orNumber;
	cout << "\nTotal scattering energy = " << D_tot/**NRM*/;
	cout << "\nTotal incoming energy = " << incomingEnergy;
	cout << "\nAveraged cross section = " << incomingEnergy*NRM;
	cout << "\nAll done. Please, press ENTER.";
}
