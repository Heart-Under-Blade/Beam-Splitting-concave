#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>
#include <iomanip>

#include "CalcTimer.h"
#include "macro.h"

#include "Mueller.hpp"
#include "CertainAggregate.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "Beam.h"
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
#include "TracerPOTotal.h"
#include "TracerBackScatterPoint.h"
#include "HandlerBackScatterPoint.h"
#include "ArgPP.h"
#include "Tracks.h"

#include "HandlerBackscatterPoint.h"
#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"
#include "HandlerPOTotal.h"

#include "ArgumentParser.h"

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

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
	parser.AddRule("random", 2, true); /* random orientarion (beta number,
gamma number)*/
	parser.AddRule("go", 0, true); // geometrical optics method
	parser.AddRule("po", 0, true); // phisical optics method
	parser.AddRule("w", 1, true); // wavelength
	parser.AddRule("b", 2, true, "po"); // beta range (begin, end)
	parser.AddRule("g", 2, true, "po"); // gamma range (begin, end)
	parser.AddRule("grid", 3, true); // backscattering grid (radius, phi, theta)
	parser.AddRule("point", zero, true, "po"); // calculate only backscatter point
	parser.AddRule("tr", 1, true); // file with trajectories
	parser.AddRule("all", 0, true); // calculate all trajectories
	parser.AddRule("abs", 1, true, "w"); // accounting of absorbtion
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
			min = Orientation::DegToRad(parser.GetDoubleValue(key, 0));
			max = Orientation::DegToRad(parser.GetDoubleValue(key, 1));
		}
		else
		{
			min = 0;
			max = particle->GetSymmetry().zenith;
		}
	}
	else if (key == "g")
	{
		number = parser.GetIntValue("random", 1);

		if (parser.IsCatched("g"))
		{
			min = Orientation::DegToRad(parser.GetDoubleValue(key, 0));
			max = Orientation::DegToRad(parser.GetDoubleValue(key, 1));
		}
		else
		{
			min = 0;
			max = particle->GetSymmetry().azimuth;
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

	bool isNew = false;
	bool isAbs = args.IsCatched("abs");

	if (isAbs)
	{
		std::string absType = args.GetStringValue("abs");

		if (absType == "new")
		{
			isAbs = true;
		}
		else if (absType == "old")
		{
			isAbs = false;
		}
		else
		{
			assert(false && "ERROR! Incorrect type of absorption.");
		}
	}

	double re = args.GetDoubleValue("ri", 0);
	double im = args.GetDoubleValue("ri", 1);
	complex refrIndex = complex(re, im);

	// TODO: AggregateBuilder

	Particle *particle = nullptr;

	std::string additionalSummary;

	additionalSummary += "Particle:\n\tType: ";

	if (args.IsCatched("pf"))
	{
		std::string filename = args.GetStringValue("p");
		particle = new Particle();
		particle->SetFromFile(filename);

		double origDMax = particle->MaximalDimension();
		cout << "from file: " << filename << endl;
		additionalSummary += "\n\nOriginal Dmax: " + std::to_string(origDMax);

		if (args.IsCatched("rs"))
		{
			double newSize = args.GetDoubleValue("rs");
			particle->Resize(newSize);
		}

		double newDMax = particle->MaximalDimension();
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
		case ParticleType::RegularColumn:
			particle = new RegularColumn(size);
			break;
		case ParticleType::Bullet:
			sup = (size.diameter*sqrt(3)*tan(Orientation::DegToRad(62)))/4;
			particle = new Bullet(size, sup);
			break;
		case ParticleType::BulletRosette:
			sup = (size.diameter*sqrt(3)*tan(Orientation::DegToRad(62)))/4;
			particle = new BulletRosette(size, sup);
			break;
//		case ParticleType::TiltedHexagonal:
//			sup = parser.argToValue<double>(vec[3]);
//			particle = new TiltedHexagonal(r, hh, ri, sup);
//			break;
		case ParticleType::ConcaveHexagonal:
			sup = args.GetDoubleValue("p", 3);
			particle = new DistortedColumn(size, sup);
			break;
		case ParticleType::HollowColumn:
			sup = args.GetDoubleValue("p", 3);
			particle = new HollowColumn(size, sup);
			break;
		case ParticleType::HexagonalAggregate:
			num = args.GetIntValue("p", 3);
			particle = new HexagonalAggregate(size, num);
			break;
		case ParticleType::CertainAggregate:
			sup = args.GetDoubleValue("p", 3);
			particle = new CertainAggregate();
			break;
		default:
			assert(false && "ERROR! Incorrect type of particle.");
			break;
		}
	}

	additionalSummary += "\n\tArea: " + to_string(particle->Area());
	additionalSummary += "\n\tNumber of facets: " + to_string(particle->nFacets);
	additionalSummary += "\n\nRefractive index: " + to_string(re);

	if (fabs(im) > FLT_EPSILON)
	{
		additionalSummary += " + i" + to_string(im);
	}

	cout << endl;
	cout << "Area: " << particle->Area() << endl;

	int maxActNo = args.GetDoubleValue("n");
	additionalSummary += "\nNumber of secondary reflections: " + to_string(reflNum);

	Scattering *scattering;
	Light incidentLight = SetIncidentLight(particle);

	if (particle->IsNonConvex())
	{
		scattering = new ScatteringNonConvex(particle, incidentLight, maxActNo,
											 refrIndex);
	}
	else
	{
		scattering = new ScatteringConvex(particle, incidentLight, maxActNo,
										  refrIndex);
	}

	Tracks trackGroups;

	string dirName = (args.IsCatched("o")) ? args.GetStringValue("o") : "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
	additionalSummary += "\nWavelength (um): " + to_string(wave);

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");

		trackGroups.ImportTracks(particle->nFacets, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}

	additionalSummary += "Method: ";

	if (args.IsCatched("po"))
	{
		additionalSummary += "Physical optics";

		Tracer *tracer;
		HandlerPO *handler;

		if (args.IsCatched("fixed"))
		{
			ScatteringSphere sphere = SetConus(args);
			additionalSummary += ", fixed orientation" << endl << endl;

			tracer = new TracerPO(particle, scattering, dirName);
			tracer->m_log = additionalSummary;
			handler = new HandlerPO(particle, &incidentLight, wave);
			handler->SetScatteringSphere(sphere);
			handler->SetTracks(&trackGroups);

			if (isAbs)
			{
				handler->EnableAbsorption(isNew);
			}

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);

			additionalSummary += "Orientation: zenith - " +
					to_string(beta) + "°, azimuth - " + to_string(gamma) +
					"°\n\n";
			tracer.m_log = additionalSummary;
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";

			double normIndex;
			OrientationRange range = args.GetRange(particle->GetSymmetry());

			if (args.IsCatched("point"))
			{
				TracerBackScatterPoint tracer(particle, scattering, dirName);
				tracer.m_log = additionalSummary;

				handler = new HandlerBackScatterPoint
						(particle, &tracer.m_incidentLight, wave);

				normIndex = range.step.azimuth /
						(range.to.azimuth - range.from.azimuth);
			}
			else
			{
				TracerPO *tracer;
				ScatteringSphere conus = SetConus(args);

				if (args.IsCatched("all"))
				{
					tracer = new TracerPOTotal(particle, scattering, dirName);
					trackGroups.push_back(TrackGroup());
					handler = new HandlerPOTotal(particle, &incidentLight, wave);
				}
				else
				{
					tracer = new TracerPO(particle, scattering, dirName);
					handler = new HandlerPO (particle, &incidentLight, wave);
				}

				normIndex = 2 * range.nAzimuth;
				handler->SetScatteringSphere(conus);
			}

			handler->SetNormIndex(normIndex);
			handler->SetTracks(&trackGroups);

			if (isAbs)
			{
				handler->EnableAbsorption(isNew);
			}

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);
			additionalSummary += "\nGrid: " + to_string(beta.number) + "x" +
					to_string(gamma.number) + "\n\n";
			tracer->m_log = additionalSummary;
			cout << additionalSummary;
			tracer->TraceRandom(range);
		}
		else
		{
			cout << endl << "error";
		}
	}
	else // go
	{
		additionalSummary += "Geometrical optics";

		TracerGO tracer(particle, scattering, dirName);
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

		if (isAbs)
		{
			handler->EnableAbsorption(isNew);
		}

		tracer.SetHandler(handler);

		if (args.IsCatched("fixed"))
		{
			Orientation orient;
			orient.zenith = args.GetDoubleValue("fixed", 0);
			orient.azimuth = args.GetDoubleValue("fixed", 1);

			additionalSummary += ", fixed orientation";
			additionalSummary += "Orientation: zenith - " +
					to_string(beta) + "°, azimuth - " + to_string(gamma) +
					"°\n\n";
			cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceFixed(orient);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";
			OrientationRange range = parser.GetRange(particle->GetSymmetry());
			additionalSummary += "\nGrid: " + to_string(beta.number) + "x" +
					to_string(gamma.number) + "\n\n";
			cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceRandom(range);
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
