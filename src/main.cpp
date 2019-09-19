#include <iostream>
#include <iomanip>
#include <fstream>
#include <chrono>
#include <assert.h>
#include <math.h>
#include <float.h>

#include "Sphere.h"
#include "RegularColumn.h"
#include "HollowColumn.h"
#include "Bullet.h"
#include "BulletRosette.h"
#include "CertainAggregate.h"
#include "DistortedColumn.h"
#include "HexagonalAggregate.h"

#include "Mueller.hpp"
#include "Beam.h"
#include "Bullet.h"
#include "common.h"
#include "PhysMtr.hpp"
#include "ArgPP.h"
#include "Tracks.h"
#include "ArgumentParser.h"

#include "HandlerTotalGO.h"
#include "HandlerTracksGO.h"
#include "HandlerPOTotal.h"
#include "HandlerBackscatterPoint.h"

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

#include "TracerGO.h"
#include "TracerPOTotal.h"
#include "TracerBackScatterPoint.h"

#include "macro.h"

#ifdef _OUTPUT_NRG_CONV
ofstream energyFile("energy.dat", ios::out);
double SSconfined=0;
int bcount=0;
#endif

using namespace std::chrono;

enum class ParticleType : int
{
	Sphere = 0,
	RegularColumn = 1,
	Bullet = 2,
	BulletRosette = 3,
	HollowColumn = 10,
	DistortedColumn = 11,
	HexagonalAggregate = 12,
	CertainAggregate = 999
};

Tracks trackGroups;

ScatteringSphere SetScatteringSphere(ArgPP &parser)
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
		std::cerr << "Error! " << __FUNCTION__;
		throw std::exception();
	}

	return AngleRange(min, max, number);
}

int main(int argc, const char* argv[])
{
	if (argc <= 1) // no arguments
	{
		std::cout << "No arguments. Press any key to exit..." << std::endl;
		getchar();
		return 1;
	}

	ArgPP args;
	ArgumentParser parser(&args);
	parser.SetArgRules();
	args.Parse(argc, argv);

	bool isNew = false;
	bool isAbs = args.IsCatched("abs");

	if (isAbs)
	{
		std::string absMethod = args.GetStringValue("abs");

		if (absMethod == "new")
		{
			isNew = true;
		}
		else if (absMethod == "old")
		{
			isNew = false;
		}
		else
		{
			assert(false && "ERROR! Incorrect method of absorption.");
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
		std::cout << "from file: " << filename << std::endl;
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

		Size size;
		size.height = args.GetDoubleValue("p", 1);
		size.diameter = args.GetDoubleValue("p", 2);

		double sup;
		int num;
		double nPar, nMer, rad;

		switch (type)
		{
		case ParticleType::Sphere:
			nPar = args.GetIntValue("p", 1);
			nMer = args.GetIntValue("p", 2);
			rad = args.GetDoubleValue("p", 3);
			particle = new Sphere(nPar, nMer, rad);
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

	additionalSummary += "\n\tArea: " + std::to_string(particle->Area());
	additionalSummary += "\n\tNumber of facets: " +
			std::to_string(particle->nElems);
	additionalSummary += "\n\nRefractive index: " + std::to_string(re);

	if (fabs(im) > FLT_EPSILON)
	{
		additionalSummary += " + i" + std::to_string(im);
	}

	int maxActNo = args.GetDoubleValue("n");
	additionalSummary += "\nNumber of secondary reflections: " + std::to_string(maxActNo);

	Scattering *scattering;

	if (particle->IsNonConvex())
	{
		scattering = new ScatteringNonConvex(refrIndex, maxActNo);
	}
	else
	{
		scattering = new ScatteringConvex(refrIndex, maxActNo);
	}

	scattering->SetActiveParticle(particle);

	Tracks trackGroups;

	std::string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
												: "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
	additionalSummary += "\nWavelength (um): " + std::to_string(wave);

	if (args.IsCatched("tr"))
	{
		std::string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}

	additionalSummary += "Method: ";

	if (args.IsCatched("po"))
	{
		additionalSummary += "Physical optics";
		TracerPO *tracer;
		HandlerPO *handler;

		if (args.IsCatched("fixed"))
		{
			additionalSummary += ", fixed orientation";
			ScatteringSphere sphere = SetScatteringSphere(args);

			if (args.IsCatched("all"))
			{
				tracer = new TracerPOTotal(particle, scattering, dirName);
				trackGroups.push_back(TrackGroup());
				handler = new HandlerPOTotal(scattering, wave);
			}
			else
			{
				tracer = new TracerPO(particle, scattering, dirName);
				handler = new HandlerPO(scattering, wave);
			}

			tracer->m_log = additionalSummary;

			handler->SetScatteringSphere(sphere);
			handler->SetTracks(&trackGroups);

			if (isAbs)
			{
				handler->EnableAbsorption(isNew);
			}

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);

			Orientation orientation;
			orientation.zenith  = args.GetDoubleValue("fixed", 0);
			orientation.azimuth = args.GetDoubleValue("fixed", 1);

			additionalSummary += "Orientation: zenith - " +
					std::to_string(orientation.zenith) + "°, azimuth - " +
					std::to_string(orientation.azimuth) + "°\n\n";
			tracer->m_log += additionalSummary;

			tracer->TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";

			double normIndex;
			OrientationRange range = parser.GetRange(particle->GetSymmetry());

			if (args.IsCatched("point"))
			{
				tracer = new TracerBackScatterPoint(particle, scattering,
													dirName);
				tracer->m_log += additionalSummary;

				handler = new HandlerBackScatterPoint(scattering, wave);

				normIndex = range.step.azimuth /
						(range.to.azimuth - range.from.azimuth);
			}
			else
			{
				ScatteringSphere conus = SetScatteringSphere(args);

				if (args.IsCatched("all"))
				{
					tracer = new TracerPOTotal(particle, scattering, dirName);
					trackGroups.push_back(TrackGroup());
					handler = new HandlerPOTotal(scattering, wave);
				}
				else
				{
					tracer = new TracerPO(particle, scattering, dirName);
					handler = new HandlerPO(scattering, wave);
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
			additionalSummary += "\nGrid: " +
					std::to_string(range.nZenith) + "x" +
					std::to_string(range.nAzimuth) + "\n\n";
			tracer->m_log = additionalSummary;
			std::cout << additionalSummary;
			tracer->TraceRandom(range);
		}
		else
		{
			std::cout << std::endl << "error";
		}
	}
	else // go
	{
		additionalSummary += "Geometrical optics";

		TracerGO tracer(particle, scattering, dirName);
		tracer.m_log = additionalSummary;
		tracer.SetIsOutputGroups(isOutputGroups);

		HandlerGO *handler;

		if (args.IsCatched("tr"))
		{
			handler = new HandlerTracksGO(scattering, wave);
			handler->SetTracks(&trackGroups);
		}
		else
		{
			handler = new HandlerTotalGO(scattering, wave);
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
			additionalSummary +=
					"Orientation: zenith - " + std::to_string(orient.zenith) +
					"°, azimuth - " + std::to_string(orient.azimuth) + "°\n\n";
			std::cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceFixed(orient);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";
			OrientationRange range = parser.GetRange(particle->GetSymmetry());
			additionalSummary += "\nGrid: " +
					std::to_string(range.nZenith) + "x" +
					std::to_string(range.nAzimuth) + "\n\n";
			std::cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceRandom(range);
		}

		delete handler;
	}

	std::cout << std::endl << "done";

	if (!args.IsCatched("close"))
	{
		getchar();
	}

	return 0;
}
