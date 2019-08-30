#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>
#include <iomanip>

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
	RegularColumn = 1,
	Bullet = 2,
	BulletRosette = 3,
	HollowColumn = 10,
	DistortedColumn = 11,
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


Light SetIncidentLight(Particle *particle)
{
	Light incidentLight;
	incidentLight.direction = Point3f(0, 0, -1);
	incidentLight.polarizationBasis = Point3f(0, 1, 0);

	Point3f point = incidentLight.direction * particle->LongRadius();
	incidentLight.direction.d_param = Point3f::DotProduct
			(point, incidentLight.direction);

	return incidentLight;
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
	ArgumentParser parser(&args);

	parser.SetArgRules();
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

	cout << "Particle: ";

	std::string additionalSummary;

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
			double ratio = newSize/particle->MaximalDimension();
			particle->Scale(ratio);
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
		case ParticleType::DistortedColumn:
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

	cout << "Refractive index: " << re;

	if (fabs(im) > FLT_EPSILON)
	{
		cout << " + i" << im;
	}

	cout << endl;
	cout << "Area: " << particle->Area() << endl;

	int maxActNo = args.GetDoubleValue("n");

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
	cout << "Wavelength (um): " << wave << endl;

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}

	cout << "Method: ";

	if (args.IsCatched("po"))
	{
		Tracer *tracer;
		HandlerPO *handler;

		if (args.IsCatched("fixed"))
		{
			ScatteringSphere sphere = SetConus(args);
			cout << ", fixed orientation" << endl << endl;

			tracer = new TracerPO(particle, scattering, dirName);
			handler = new HandlerPO(particle, &incidentLight, wave);
			handler->SetScatteringSphere(sphere);
			handler->SetTracks(&trackGroups);

			if (isAbs)
			{
				handler->EnableAbsorption(isNew);
			}

			tracer->m_summary = additionalSummary;

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);

			Orientation orientation;
			orientation.zenith  = args.GetDoubleValue("fixed", 0);
			orientation.azimuth = args.GetDoubleValue("fixed", 1);
			tracer->TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
			cout << ", random orientation" << endl << endl;

			double normIndex;
			OrientationRange range = parser.GetRange(particle->GetSymmetry());

			HandlerPO *handler;

			if (args.IsCatched("point"))
			{
				TracerBackScatterPoint tracer(particle, scattering, dirName);
				tracer.m_summary = additionalSummary;

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

				tracer->m_summary = additionalSummary;

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
			tracer->TraceRandom(range);
		}
		else
		{
			cout << endl << "error";
		}

		delete handler;
		delete tracer;
	}
	else // GO
	{
		cout << "Geometrical optics";

		TracerGO tracer(particle, scattering, dirName);
		tracer.m_summary = additionalSummary;
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

		if (isAbs)
		{
			handler->EnableAbsorption(isNew);
		}

		tracer.SetHandler(handler);

		if (args.IsCatched("fixed"))
		{
			cout << ", fixed orientation, ";
			Orientation orient;
			orient.zenith = args.GetDoubleValue("fixed", 0);
			orient.azimuth = args.GetDoubleValue("fixed", 1);
			cout << "zenith " << orient.zenith
				 << "°, azimuth " << orient.azimuth << "°" << endl << endl;
			tracer.TraceFixed(orient);
		}
		else if (args.IsCatched("random"))
		{
			cout << ", random orientation, ";
			OrientationRange range = parser.GetRange(particle->GetSymmetry());
			cout << "grid: " << range.nZenith << "x" << range.nAzimuth
				 << endl << endl;
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
