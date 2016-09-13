#include <iostream>
#include <time.h>

//#include "test.h"

#include <float.h>
#include "global.h"

#include "Mueller.hpp"
#include "Hexagon.h"
#include "Tracing.h"
#include "Beam.h"
#include "PhysMtr.hpp"

matrix back(4,4),	///< Mueller matrix in backward direction
		forw(4,4);	///< Mueller matrix in forward direction

Arr2D mxd(0, 0, 0, 0);
double SizeBin;
double gammaNorm, betaNorm;
double betaDistrProbability;
double incomingEnergy;
clock_t totalTime;
Point3f incidentDir(0, 0, 1);
Point3f polarizationBasis(0, 1, 0);

long long ddddd = 0;

void HandleBeams(std::vector<OutBeam> &outBeams)
{
	for (unsigned int i = 0; i < outBeams.size(); ++i)
	{
		++ddddd;
		Beam &bm = outBeams.at(i).beam;
		bm.RotateSpherical(incidentDir, polarizationBasis);

		unsigned int szP = outBeams.at(i).trackSize;
		// If there is a special list of the beams to take into account, check if the beam "bm" is on the list
/// TODO: для отдельных траекторий
//		if (isSorting)
//		{
//			bool flag = false;
//			std::list<Chain>::const_iterator c = mask.begin();

//			for (; c!= mask.end(); ++c)
//			{
//				std::list<unsigned int >::const_iterator it = c->Ch.begin();

//				if (szP!=c->sz)
//				{
//					continue;
//				}

//				std::list<unsigned int >::const_iterator fs = bm.BeginP();

//				for (;it!=c->Ch.end() && (*it)==(*fs); ++it, ++fs);

//				if (it==c->Ch.end())
//				{
//					flag = true;
//					break;
//				}
//			}

//			if (!flag)
//			{
//				return;
//			}
//		}

		double cross = bm.CrossSection();
		double Area = betaDistrProbability*cross;
		matrix bf = Mueller(bm.JMatrix);
		//----------------------------------------------------------------------------
		// Collect the beam in array

		if(bm.direction.Z >= 1-DBL_EPSILON)
		{
			back += Area*bf;
		}
		else
		{
			if(bm.direction.Z <= DBL_EPSILON-1)
			{
				forw += Area*bf;
			}
			else
			{
				// Rotate the Mueller matrix of the beam to appropriate coordinate system
				const unsigned int ZenAng = round(acos(bm.direction.Z)/SizeBin);

				double tmp = SQR(bm.direction.Y);

				if (tmp > DBL_EPSILON)
				{
					tmp = acos(bm.direction.X/sqrt(SQR(bm.direction.X)+tmp));

					if (bm.direction.Y < 0)
					{
						tmp = M_2PI-tmp;
					}

					tmp *= -2.0;
					double cs = cos(tmp),
							sn = sin(tmp);

					RightRotateMueller(bf, cs, sn);
				}

				mxd.insert(0, ZenAng, Area*bf);
			}
		}
	}
}

void Trace(double beta, double gamma, Tracing &tracer)
{
	if (ddddd == 579)
	{
		int ff = 0;
	}
	double square = tracer.TraceParticle(beta, gamma); /// TODO: трассировка
	incomingEnergy += betaDistrProbability * square;
//	tracer.ClearTracing();
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

void WriteStatisticsToConsole(int orNumber, double D_tot, double NRM)
{
	using namespace std;
	cout << endl << ddddd << endl;
	cout << "\nTotal number of body orientation = " << orNumber;
	cout << "\nTotal scattering energy = " << D_tot/**NRM*/;
	cout << "\nTotal incoming energy = " << incomingEnergy;
	cout << "\nAveraged cross section = " << incomingEnergy*NRM;
	cout << "\nAll done. Please, press ENTER.";
	getchar();
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
			frw << "\n" << f[0]*NRM << " " << f[1]/f[0] << " " << f[1]/f[0] << " " << f[2]/f[0];
		}

		if (b[0] <= DBL_EPSILON)
		{
			bck << "\n0 0 0 0";
		}
		else
		{
			bck << "\n" << b[0]*NRM << " " << b[1]/b[0] << " " << -b[1]/b[0] << " " << b[2]/b[0];
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

void Calculate()
{
	/// params
	int particleType = 0;
	double radius = 40;
	double halfHeight = 100;
	complex refractionIndex = complex(1.31, 0.0);
	int orNumber_gamma = 101;
	int orNumber_beta = 100;
	int ThetaNumber = 180;
	int interReflNum = 0;
	bool isRandom = false;

	bool isOpticalPath = true;
	bool isSorting = false;
	int EDF = 0;

	Particle *particle = nullptr;
	switch (particleType) {
	case 0:
		particle = new Hexagon(radius, halfHeight, refractionIndex);
		betaNorm = M_PI/(2.0*orNumber_beta);
		break;
	case 1:
		/// TODO ...
		break;
	default:
		break;
	}
//	CreateParticle(particle, particleType,
//				   halfHeight, radius, orNumber_beta, refractionIndex);

	gammaNorm = M_PI/(3.0*orNumber_gamma);

	double beta, gamma;
	incomingEnergy = 0;
	betaDistrProbability = 0;

	srand(time(NULL));

	// the arrays for exact backscattering and forwardscattering Mueller matrices
	back.Fill(0);
	forw.Fill(0);

	SizeBin = M_PI/ThetaNumber; // the size of the bin (radians)

	mxd = Arr2D(1, ThetaNumber + 1, 4, 4);
	mxd.ClearArr();

	Tracing tracer(particle, isOpticalPath, incidentDir,
				   polarizationBasis, interReflNum);

	clock_t timer = clock();

	if (isRandom)
	{
		gamma = 2.0*M_PI * ((double)(rand()) / ((double)RAND_MAX + 1.0));
		beta = acos(((double)(rand()) / ((double)RAND_MAX))*2.0 - 1.0);
		betaDistrProbability = 1.0;

		for (int i = 0; i < orNumber_beta; ++i)
		{
			std::cout << i << "/" << orNumber_beta << std::endl;

			for (int j = 0; j < orNumber_gamma; ++j)
			{
				Trace(beta, gamma, tracer);
			}
		}
	}
	else
	{
//		beta = 20;
//		gamma = 0;
//		Trace(beta, gamma, tracer);

		for (int i = 0; i < orNumber_beta; ++i)
		{
			std::cout << i << "/" << orNumber_beta << std::endl;

			beta = (i + 0.5)*betaNorm;
			betaDistrProbability = sin(beta);

			for (int j = 0; j < orNumber_gamma; ++j)
			{
//				if (ddddd == 583)
//				{
//					int ff = 0;
//				}
				gamma = (j + 0.5)*gammaNorm;
				Trace(beta, gamma, tracer);
			}
		}
	}

	HandleBeams(tracer.outcomingBeams);

	timer = clock() - timer;
	totalTime = timer/CLOCKS_PER_SEC;

	std::cout << "\nTotal time of calculation = " << totalTime << " seconds";

	//Integrating
//	double aa = back[0][0];
//	double bb = forw[0][0];
	double D_tot = back[0][0] + forw[0][0];

	for (int j = 0; j <= ThetaNumber; ++j)
	{
//		std::cout << j << ": " << D_tot << ", " << mxd(0, j, 0, 0) << std::endl;
		D_tot += mxd(0, j, 0, 0);
	}

	//Normalizing coefficient
	double orNumber = orNumber_gamma*orNumber_beta;
	double NRM;

	if (!isRandom)
	{
		NRM = M_PI/((double)orNumber*2.0);
	}
	else
	{
		NRM = 1.0/(double)orNumber;
	}
	//const double NRM = 1.0/(4.0*M_PI);

	ExtractPeaks(EDF, NRM, ThetaNumber);

	WriteResultsToFile(ThetaNumber, NRM);

	WriteStatisticsToFile(timer, orNumber, D_tot, NRM);

	WriteStatisticsToConsole(orNumber, D_tot, NRM);
	getchar();

	delete particle;
}

int main()
{
//	testHexagonBuilding();
//	testHexagonRotate();
	Calculate();
	return 0;
}

