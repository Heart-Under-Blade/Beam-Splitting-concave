#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>

#include "CalcTimer.h"
#include "macro.h"

#include "Mueller.hpp"
#include "RegularColumn.h"
#include "HollowColumn.h"
#include "HexagonalAggregate.h"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "DistortedColumn.h"
#include "common.h"
#include "PhysMtr.hpp"
#include "TracerGO.h"
#include "TracerBackScatterPoint.h"
#include "ArgPP.h"
#include "Tracks.h"
#include "Handler.h"

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std;
using namespace chrono;

enum class ParticleType : int
{
	RegularColumn = 1,
	Bullet = 2,
	BulletRosette = 3,
	HollowColumn = 10,
	DistortedColumn = 11,
	HexagonalAggregate = 12,
	CertainAggregate = 999
};

void SetArgRules(ArgPP &parser)
{
	int zero = 0;
	parser.AddRule("p", '+'); // particle (type, size, ...)
	parser.AddRule("ri", 2); // refractive index (Re and Im parts)
	parser.AddRule("n", 1); // number of internal reflection
	parser.AddRule("pf", zero, true); // particle (filename)
	parser.AddRule("fixed", 2, true); // fixed orientarion (beta, gamma)
	parser.AddRule("random", 2, true); // random orientarion (beta number, gamma number)
	parser.AddRule("go", 0, true); // geometrical optics method
	parser.AddRule("po", 0, true); // phisical optics method
	parser.AddRule("w", 1, true); // wavelength
	parser.AddRule("b", 2, true); // beta range (begin, end)
	parser.AddRule("g", 2, true); // gamma range (begin, end)
	parser.AddRule("conus", 3, true, "po"); // calculate only backscatter cone (radius, phi, theta)
	parser.AddRule("point", zero, true, "po"); // calculate only backscatter poin
	parser.AddRule("con20", zero, true, "po");
	parser.AddRule("gr", zero, true);
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

		if (parser.IsCatched(key))
		{
			min = Angle3d::DegToRad(parser.GetDoubleValue(key, 0));
			max = Angle3d::DegToRad(parser.GetDoubleValue(key, 1));
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

		if (parser.IsCatched(key))
		{
			min = Angle3d::DegToRad(parser.GetDoubleValue(key, 0));
			max = Angle3d::DegToRad(parser.GetDoubleValue(key, 1));
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

Light SetIncidentLight(Particle *particle)
{
	Light incidentLight;
	incidentLight.direction = Point3f(0, 0, -1);
	incidentLight.polarizationBasis = Point3f(0, 1, 0);

	Point3f point = incidentLight.direction * particle->ComputeRotationRadius();
	incidentLight.direction.d_param = Point3f::DotProduct(point, incidentLight.direction);
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

	if (args.IsCatched("pf"))
	{
		std::string filename = args.GetStringValue("p");
		particle = new Particle();
		particle->SetFromFile(filename);
		particle->SetRefractiveIndex(complex(refrIndex));
	}
	else
	{
		ParticleType type = (ParticleType)args.GetIntValue("p", 0);
		Size size;
		size.height = args.GetDoubleValue("p", 1);
		size.diameter = args.GetDoubleValue("p", 2);

		double sup;
		int num;

		switch (type)
		{
		case ParticleType::RegularColumn:
			particle = new RegularColumn(refrIndex, size);
			break;
		case ParticleType::Bullet:
			sup = (size.diameter*sqrt(3)*tan(Angle3d::DegToRad(62)))/4;
			particle = new Bullet(refrIndex, size, sup);
			break;
		case ParticleType::BulletRosette:
			sup = (size.diameter*sqrt(3)*tan(Angle3d::DegToRad(62)))/4;
			particle = new BulletRosette(refrIndex, size, sup);
			break;
		case ParticleType::DistortedColumn:
			sup = args.GetDoubleValue("p", 3);
			particle = new DistortedColumn(refrIndex, size, sup);
			break;
		case ParticleType::HollowColumn:
			sup = args.GetDoubleValue("p", 3);
			particle = new HollowColumn(refrIndex, size, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = args.GetIntValue("p", 3);
			particle = new HexagonalAggregate(refrIndex, size, num);
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

	int maxActNo = args.GetDoubleValue("n");

	Scattering *scattering;
	Light incidentLight = SetIncidentLight(particle);

	if (particle->IsNonConvex())
	{
		scattering = new ScatteringNonConvex(particle, incidentLight, maxActNo);
	}
	else
	{
		scattering = new ScatteringConvex(particle, incidentLight, maxActNo);
	}

	Tracks trackGroups;

	string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
										   : "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w")
									  : 0;

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle->nElems, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}

	if (args.IsCatched("po"))
	{
		if (args.IsCatched("fixed"))
		{
			Conus bsCone = SetCone(args);
			TracerPO tracer(particle, scattering, dirName);

			HandlerPO *handler = new HandlerPO(particle, &incidentLight, wave);
			handler->SetTracks(&trackGroups);
			handler->SetScatteringConus(bsCone);
			handler->SetAbsorbtionAccounting(isAbs);

			tracer.SetIsOutputGroups(isOutputGroups);
			tracer.SetHandler(handler);

			Angle3d orientation;
			orientation.beta  = args.GetDoubleValue("fixed", 0);
			orientation.gamma = args.GetDoubleValue("fixed", 1);
			tracer.TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
			HandlerPO *handler;
			LightTracer *tracer;

			if (args.IsCatched("point"))
			{
				tracer = new TracerBackScatterPoint(particle, scattering, dirName);
				handler = new HandlerBackScatterPoint(particle, &incidentLight, wave);

				if (args.IsCatched("con20"))
				{
					handler->setCon20(false);
				}
			}
			else
			{
				tracer = new TracerPO(particle, scattering, dirName);
				handler = new HandlerPO(particle, &incidentLight, wave);

				Conus bsCone = SetCone(args);
				handler->SetScatteringConus(bsCone);
			}

			AngleRange zenith = GetRange(args, "b", particle);
			AngleRange azimuth = GetRange(args, "g", particle);

			double normIndex = azimuth.step/azimuth.norm;
			handler->SetNormIndex(normIndex);
			handler->SetTracks(&trackGroups);
			handler->SetAbsorbtionAccounting(isAbs);

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);
			tracer->TraceRandom(zenith, azimuth);

			delete handler;
			delete tracer;
		}
		else
		{
			cout << endl << "error";
		}
	}
	else // GO
	{
		TracerGO tracer(particle, scattering, dirName);
		tracer.SetIsOutputGroups(isOutputGroups);

		HandlerGO *handler;

		if (args.IsCatched("tr"))
		{
			handler = new HandlerTracksGO(particle, &incidentLight, wave);
			handler->SetTracks(&trackGroups);
		}
		else
		{
			handler = new HandlerTotalGO(particle, &incidentLight, wave);
		}

		handler->SetAbsorbtionAccounting(isAbs);
		tracer.SetHandler(handler);

		if (args.IsCatched("fixed"))
		{
			Angle3d orientation;
			orientation.beta  = args.GetDoubleValue("fixed", 0);
			orientation.gamma = args.GetDoubleValue("fixed", 1);
			tracer.TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
			AngleRange zenith = GetRange(args, "b", particle);
			AngleRange azimuth = GetRange(args, "g", particle);
			tracer.TraceRandom(zenith, azimuth);
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
