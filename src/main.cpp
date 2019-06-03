#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>
#include <iomanip>

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
#include "TracerPOTotal.h"
#include "TracerBackScatterPoint.h"
#include "HandlerBackScatterPoint.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "HandlerPOTotal.h"
#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"

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
	parser.AddRule("pf", zero, true); // particle (filename)
	parser.AddRule("rs", 1, true, "pf"); // resize particle (new size)
	parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
	parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
	parser.AddRule("go", 0, true); // geometrical optics method
	parser.AddRule("po", 0, true); // phisical optics method
	parser.AddRule("w", 1, true); // wavelength
	parser.AddRule("b", 2, true, "po"); // beta range (begin, end)
	parser.AddRule("g", 2, true, "po"); // gamma range (begin, end)
	parser.AddRule("grid", 3, true); // backscattering grid (radius, phi, theta)
	parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
	parser.AddRule("tr", 1, true); // file with trajectories
	parser.AddRule("all", 0, true); // calculate all trajectories
	parser.AddRule("abs", zero, true, "w"); // accounting of absorbtion
	parser.AddRule("close", 0, true); // closing of program after calculation
	parser.AddRule("o", 1, true); // output folder name
}

ScatteringSphere SetConus(ArgPP &parser)
{
	double radius = parser.GetDoubleValue("grid", 0);
	int nPhi = parser.GetDoubleValue("grid", 1);
	int nTheta = parser.GetDoubleValue("grid", 2);
	return ScatteringSphere(radius, nPhi, nTheta);
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
		cout << "No arguments. Press any key to exit..." << endl;
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

	cout << "Particle: ";

	std::string additionalSummary;

	if (args.IsCatched("pf"))
	{
		std::string filename = args.GetStringValue("p");
		particle = new Particle();
		particle->SetFromFile(filename);
		particle->SetRefractiveIndex(complex(refrIndex));

		double origDMax = particle->MaximalDimention();
		cout << "from file: " << filename << endl;
		additionalSummary += "\n\nOriginal Dmax: " + std::to_string(origDMax);

		if (args.IsCatched("rs"))
		{
			double newSize = args.GetDoubleValue("rs");
			particle->Resize(newSize);
		}

		double newDMax = particle->MaximalDimention();
		additionalSummary += ", new Dmax: " + std::to_string(newDMax)
				+ ", resize factor: " + std::to_string(newDMax/origDMax) + '\n';
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

	cout << "Refractive index: " << re;

	if (fabs(im) > FLT_EPSILON)
	{
		cout << " + i" << im;
	}

	cout << endl;

//	particle->Output();
	cout << "Area:" << particle->Area() << endl;

	int reflNum = args.GetDoubleValue("n");
	cout << "Number of secondary reflections: " << reflNum << endl;

	string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
										   : "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w")
									  : 0;
	cout << "Wavelength (um): " << wave << endl;

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle->nFacets, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}


	cout << "Method: ";

	if (args.IsCatched("po"))
	{
		cout << "Physical optics";

		if (args.IsCatched("fixed"))
		{
			cout << ", fixed orientation" << endl << endl;

			ScatteringSphere bsCone = SetConus(args);

			double beta  = args.GetDoubleValue("fixed", 0);
			double gamma = args.GetDoubleValue("fixed", 1);

			TracerPO tracer(particle, reflNum, dirName);
			tracer.m_summary = additionalSummary;

			HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight, wave);
			handler->SetTracks(&trackGroups);
			handler->SetScatteringSphere(bsCone);
			handler->SetAbsorptionAccounting(isAbs);

			tracer.SetIsOutputGroups(isOutputGroups);
			tracer.SetHandler(handler);
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			cout << ", random orientation" << endl << endl;
			AngleRange beta = GetRange(args, "b", particle);
			AngleRange gamma = GetRange(args, "g", particle);

			HandlerPO *handler;

			if (args.IsCatched("point"))
			{
				TracerBackScatterPoint tracer(particle, reflNum, dirName);
				tracer.m_summary = additionalSummary;

				handler = new HandlerBackScatterPoint(particle, &tracer.m_incidentLight, wave);
				double normIndex = gamma.step/gamma.norm;
				handler->SetNormIndex(normIndex);
				handler->SetTracks(&trackGroups);
				handler->SetAbsorptionAccounting(isAbs);

				tracer.SetIsOutputGroups(isOutputGroups);
				tracer.SetHandler(handler);
				tracer.TraceRandom(beta, gamma);
			}
			else
			{
				TracerPO *tracer;
				ScatteringSphere conus = SetConus(args);

				if (args.IsCatched("all"))
				{
					tracer = new TracerPOTotal(particle, reflNum, dirName);
					trackGroups.push_back(TrackGroup());
					handler = new HandlerPOTotal(particle, &tracer->m_incidentLight, wave);
				}
				else
				{
					tracer = new TracerPO(particle, reflNum, dirName);
					handler = new HandlerPO(particle, &tracer->m_incidentLight, wave);
				}

				tracer->m_summary = additionalSummary;

				handler->SetTracks(&trackGroups);
				handler->SetScatteringSphere(conus);
				handler->SetAbsorptionAccounting(isAbs);

				tracer->SetIsOutputGroups(isOutputGroups);
				tracer->SetHandler(handler);
				tracer->TraceRandom(beta, gamma);
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
		cout << "Geometrical optics";

		TracerGO tracer(particle, reflNum, dirName);
		tracer.m_summary = additionalSummary;
		tracer.SetIsOutputGroups(isOutputGroups);

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

		handler->SetAbsorptionAccounting(isAbs);
		tracer.SetHandler(handler);

		if (args.IsCatched("fixed"))
		{
			cout << ", fixed orientation, ";
			double beta  = args.GetDoubleValue("fixed", 0);
			double gamma = args.GetDoubleValue("fixed", 1);
			cout << "zenith " << beta << "°, azimuth " << gamma << "°" << endl << endl;
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			cout << ", random orientation, ";
			AngleRange beta = GetRange(args, "b", particle);
			AngleRange gamma = GetRange(args, "g", particle);
			cout << "grid: " << beta.number << "x" << gamma.number << endl << endl;
			tracer.TraceRandom(beta, gamma);
		}

		delete handler;
	}

	cout << endl << "done";

	if (!args.IsCatched("close"))
	{
		getchar();
	}

	return 0;
}
