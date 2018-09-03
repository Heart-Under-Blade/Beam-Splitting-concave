#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>

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
#include "TracerGO.h"
#include "TracerBackScatterPoint.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "Handler.h"

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

void SetArgRules(ArgPP &parser)
{
	int zero = 0;
	parser.AddRule("p", '+'); // particle (type, size, ...)
	parser.AddRule("ri", 2); // refractive index (Re and Im parts)
	parser.AddRule("n", 1); // number of internal reflection
	parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
	parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
	parser.AddRule("go", 0, true); // geometrical optics method
	parser.AddRule("po", 0, true); // phisical optics method
	parser.AddRule("w", 1, true); // wavelength
	parser.AddRule("b", 2, true, "po"); // beta range (begin, end)
	parser.AddRule("g", 2, true, "po"); // gamma range (begin, end)
	parser.AddRule("conus", 3, true, "po"); // calculate only backscatter cone (radius, phi, theta)
	parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
	parser.AddRule("tr", 1, true); // file with trajectories
	parser.AddRule("all", 0, true); // calculate all trajectories
	parser.AddRule("abs", zero, true, "w"); // accounting of absorbtion
	parser.AddRule("close", 0, true); // closing of program after calculation
	parser.AddRule("o", 1, true); // output folder name
}

Conus SetCone(ArgPP &parser)
{
	double radius = parser.GetDoubleValue("conus", 0);
	int phiCount = parser.GetDoubleValue("conus", 1);
	int thetaCount = parser.GetDoubleValue("conus", 2);
	return Conus(radius, phiCount, thetaCount);
}

AngleRange GetRange(const ArgPP &parser, const std::string &key,
					Particle *particle)
{
	int number;
	double min, max;

	if (key == "b")
	{
		number = parser.GetIntValue("random", 0);

		if (parser.IsCatched("b"))
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

		if (parser.IsCatched("g"))
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
	if (argc <= 1) // no arguments
	{
		cout << endl << "No arguments. Press any key to exit...";
		getchar();
		return 1;
	}

	ArgPP args;
	SetArgRules(args);
	args.Parse(argc, argv);

	bool isAbs = args.IsCatched("abs");

	double re = args.GetDoubleValue("ri", 0);
	double im = args.GetDoubleValue("ri", 1);
	complex refrIndex = complex(re, im);

	// TODO: AggregateBuilder

	Particle *particle = nullptr;

	if (args.GetArgNumber("p") == 1)
	{
		std::string filename = args.GetStringValue("pf");
		particle = new Particle();
		particle->SetFromFile(filename);
		particle->SetRefractiveIndex(complex(refrIndex));
	}
	else
	{
		ParticleType type = (ParticleType)args.GetIntValue("p", 0);
		double height = args.GetDoubleValue("p", 1);
		double diameter = args.GetDoubleValue("p", 2);

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
			sup = args.GetDoubleValue("p", 3);
			particle = new ConcaveHexagonal(refrIndex, diameter, height, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = args.GetIntValue("p", 3);
			particle = new HexagonalAggregate(refrIndex, diameter, height, num);
			break;
		case ParticleType::CertainAggregate:
			sup = args.GetDoubleValue("p", 3);
			particle = new CertainAggregate(refrIndex, sup);
			break;
		default:
			assert(false && "ERROR! Incorrect type of particle.");
			break;
		}
	}

	particle->Output();

	int reflNum = args.GetDoubleValue("n");
	string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
												: "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w")
									  : 0;

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle->nFacets, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}

	if (args.IsCatched("po"))
	{
		if (args.IsCatched("fixed"))
		{
			Conus bsCone = SetCone(args);

			double beta  = args.GetDoubleValue("fixed", 0);
			double gamma = args.GetDoubleValue("fixed", 1);

			TracerPO tracer(particle, reflNum, dirName);

			HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight, wave);
			handler->SetTracks(&trackGroups);
			handler->SetScatteringConus(bsCone);
			handler->SetAbsorbtionAccounting(isAbs);

			tracer.SetIsOutputGroups(isOutputGroups);
			tracer.SetHandler(handler);
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			AngleRange beta = GetRange(args, "b", particle);
			AngleRange gamma = GetRange(args, "g", particle);

			HandlerPO *handler;

			if (args.IsCatched("point"))
			{
				TracerBackScatterPoint tracer(particle, reflNum, dirName);

				handler = new HandlerBackScatterPoint(particle, &tracer.m_incidentLight, wave);
				double normIndex = gamma.step/gamma.norm;
				handler->SetNormIndex(normIndex);
				handler->SetTracks(&trackGroups);
				handler->SetAbsorbtionAccounting(isAbs);

				tracer.SetIsOutputGroups(isOutputGroups);
				tracer.SetHandler(handler);
				tracer.Trace(beta, gamma, trackGroups, wave);
			}
			else
			{
				Conus bsCone = SetCone(args);

				TracerPO tracer(particle, reflNum, dirName);

				handler = new HandlerPO(particle, &tracer.m_incidentLight, wave);
				handler->SetTracks(&trackGroups);
				handler->SetScatteringConus(bsCone);
				handler->SetAbsorbtionAccounting(isAbs);

				tracer.SetIsOutputGroups(isOutputGroups);
				tracer.SetHandler(handler);
				tracer.TraceRandom(beta, gamma);
//				tracer.TraceIntervalPO2(betaR, gammaR, bsCone, trackGroups, wave);
			}

			delete handler;
		}
		else
		{
			cout << endl << "error";
		}
	}
	else // go
	{
		TracerGO tracer(particle, reflNum, dirName);
		tracer.SetIsOutputGroups(isOutputGroups);

		AngleRange beta = GetRange(args, "b", particle);
		AngleRange gamma = GetRange(args, "g", particle);

		HandlerGO *handler;

		if (args.IsCatched("tr"))
		{
			handler = new HandlerTracksGO(particle, &tracer.m_incidentLight, wave);
			handler->SetTracks(&trackGroups);
		}
		else
		{
			handler = new HandlerTotalGO(particle, &tracer.m_incidentLight, wave);
		}

		handler->SetAbsorbtionAccounting(isAbs);
		tracer.SetHandler(handler);
		tracer.TraceRandom(beta, gamma);

		delete handler;
	}

	cout << endl << "done";

	if (!args.IsCatched("close"))
	{
		getchar();
	}

	return 0;
}
