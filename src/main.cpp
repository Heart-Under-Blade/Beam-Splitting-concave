#include <iostream>
#include <fstream>
#include <assert.h>
#include <float.h>
#include <chrono>
#include <iomanip>

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

#include "ArgumentParser.h"

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

Conus SetCone(ArgPP &parser)
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

	Point3f point = incidentLight.direction * particle->RotationRadius();
	incidentLight.direction.d_param = Point3f::DotProduct(point, incidentLight.direction);
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

	bool isAbs = args.IsCatched("abs");

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
		particle->SetRefractiveIndex(complex(refrIndex));

		double origDMax = particle->MaximalDimention();
		additionalSummary += "\tfile \"" + filename + "\"";
		additionalSummary += "\n\tOriginal Dmax: " + std::to_string(origDMax);

		if (args.IsCatched("rs"))
		{
			double newSize = args.GetDoubleValue("rs");
			double ratio = newSize/particle->MaximalDimension();
			particle->Scale(ratio);
		}

		double newDMax = particle->MaximalDimention();
		additionalSummary += "\n\tNew Dmax: " + std::to_string(newDMax)
				+ " (resize factor: " + std::to_string(newDMax/origDMax) + ')';
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
			sup = (size.diameter*sqrt(3)*tan(Orientation::DegToRad(62)))/4;
			particle = new Bullet(refrIndex, size, sup);
			break;
		case ParticleType::BulletRosette:
			sup = (size.diameter*sqrt(3)*tan(Orientation::DegToRad(62)))/4;
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
			particle = new CertainAggregate(refrIndex);
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

<<<<<<< HEAD
	int reflNum = args.GetDoubleValue("n");
	additionalSummary += "\nNumber of secondary reflections: " + to_string(reflNum);
=======
	cout << endl;

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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

	string dirName = (args.IsCatched("o")) ? args.GetStringValue("o")
										   : "M";
	bool isOutputGroups = args.IsCatched("gr");
	double wave = args.IsCatched("w") ? args.GetDoubleValue("w") : 0;
	additionalSummary += "\nWavelength (um): " + to_string(wave);

	if (args.IsCatched("tr"))
	{
		string trackFileName = args.GetStringValue("tr");
		trackGroups.ImportTracks(particle->nElems, trackFileName);
		trackGroups.shouldComputeTracksOnly = !args.IsCatched("all");
	}


	additionalSummary += "\nMethod: ";

	if (args.IsCatched("po"))
	{
		additionalSummary += "Physical optics";

		HandlerPO *handler;
		LightTracer *tracer;

		if (args.IsCatched("fixed"))
		{
			additionalSummary += ", fixed orientation";

<<<<<<< HEAD
			ScatteringSphere sphere = SetConus(args);

			double beta  = args.GetDoubleValue("fixed", 0);
			double gamma = args.GetDoubleValue("fixed", 1);

			TracerPO tracer(particle, reflNum, dirName);
			tracer.m_log = additionalSummary;

			HandlerPO *handler = new HandlerPO(particle, &tracer.m_incidentLight, wave);
			handler->SetTracks(&trackGroups);
			handler->SetScatteringSphere(sphere);
			handler->SetAbsorptionAccounting(isAbs);

			tracer.SetIsOutputGroups(isOutputGroups);
			tracer.SetHandler(handler);
			additionalSummary += "Orientation: zenith - " +
					to_string(beta) + "°, azimuth - " + to_string(gamma) +
					"°\n\n";
			tracer.m_log = additionalSummary;
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";
			AngleRange beta = GetRange(args, "b", particle);
			AngleRange gamma = GetRange(args, "g", particle);

			HandlerPO *handler;

			if (args.IsCatched("point"))
			{
				TracerBackScatterPoint tracer(particle, reflNum, dirName);
				tracer.m_log = additionalSummary;

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

				handler->SetTracks(&trackGroups);
				handler->SetScatteringSphere(conus);
				handler->SetAbsorptionAccounting(isAbs);

				tracer->SetIsOutputGroups(isOutputGroups);
				tracer->SetHandler(handler);
				additionalSummary += "\nGrid: " + to_string(beta.number) + "x" +
						to_string(gamma.number) + "\n\n";
				tracer->m_log = additionalSummary;
				cout << additionalSummary;
				tracer->TraceRandom(beta, gamma);
=======
			Conus bsCone = SetCone(args);
			tracer = new TracerPO(particle, scattering, dirName);

			handler = new HandlerPO(particle, &incidentLight, wave);
			handler->SetScatteringConus(bsCone);
			handler->SetTracks(&trackGroups);
			handler->SetAbsorbtionAccounting(isAbs);

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);

			Orientation orientation;
			orientation.zenith  = args.GetDoubleValue("fixed", 0);
			orientation.azimuth = args.GetDoubleValue("fixed", 1);
			tracer->TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
			}

			AngleRange zenith = parser.GetRange("b", particle->GetSymmetry());
			AngleRange azimuth = parser.GetRange("g", particle->GetSymmetry());

			double normIndex = azimuth.step/azimuth.norm;
			handler->SetNormIndex(normIndex);
			handler->SetTracks(&trackGroups);
			handler->SetAbsorbtionAccounting(isAbs);

			tracer->SetIsOutputGroups(isOutputGroups);
			tracer->SetHandler(handler);
			tracer->TraceRandom(zenith, azimuth);
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
<<<<<<< HEAD
		additionalSummary += "Geometrical optics";

		TracerGO tracer(particle, reflNum, dirName);
=======
		TracerGO tracer(particle, scattering, dirName);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
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

		handler->SetAbsorptionAccounting(isAbs);
		tracer.SetHandler(handler);

		if (args.IsCatched("fixed"))
		{
<<<<<<< HEAD
			additionalSummary += ", fixed orientation";
			double beta  = args.GetDoubleValue("fixed", 0);
			double gamma = args.GetDoubleValue("fixed", 1);
			additionalSummary += "Orientation: zenith - " +
					to_string(beta) + "°, azimuth - " + to_string(gamma) +
					"°\n\n";
			cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceFixed(beta, gamma);
		}
		else if (args.IsCatched("random"))
		{
			additionalSummary += ", random orientation";
			AngleRange beta = GetRange(args, "b", particle);
			AngleRange gamma = GetRange(args, "g", particle);
			additionalSummary += "\nGrid: " + to_string(beta.number) + "x" +
					to_string(gamma.number) + "\n\n";
			cout << additionalSummary;
			tracer.m_log = additionalSummary;
			tracer.TraceRandom(beta, gamma);
=======
			Orientation orientation;
			orientation.zenith  = args.GetDoubleValue("fixed", 0);
			orientation.azimuth = args.GetDoubleValue("fixed", 1);
			tracer.TraceFixed(orientation);
		}
		else if (args.IsCatched("random"))
		{
			AngleRange zenith = parser.GetRange("b", particle->GetSymmetry());
			AngleRange azimuth = parser.GetRange("g", particle->GetSymmetry());
			tracer.TraceRandom(zenith, azimuth);
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
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
