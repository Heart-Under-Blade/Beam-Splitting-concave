#include "Tracer.h"

#include <ostream>
#include <iostream>
#include <assert.h>
#include "global.h"
#include "macro.h"
//#ifdef _TRACK_ALLOW
//std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
//#endif

#include "ScatteringConvex.h"
#include "ScatteringNonConvex.h"

using namespace std;

Tracer::Tracer(Particle *particle, int maxActNo, const string &resultFileName)
	: m_particle(particle),
	  m_resultDirName(resultFileName)
{
	SetIncidentLight(particle);

	if (particle->IsNonConvex())
	{
		m_scattering = new ScatteringNonConvex(particle, m_incidentLight, maxActNo);
	}
	else
	{
		m_scattering = new ScatteringConvex(particle, m_incidentLight, maxActNo);
	}
}

Tracer::~Tracer()
{
}

void Tracer::SetIncidentLight(Particle *particle)
{
	m_incidentLight.direction = Point3f(0, 0, -1);
	m_incidentLight.polarizationBasis = Point3f(0, 1, 0);

	Point3f point = m_incidentLight.direction * particle->ComputeRotationRadius();
	m_incidentLight.direction.d_param = Point3f::DotProduct(point, m_incidentLight.direction);
}

void Tracer::OutputOrientationToLog(int i, int j, ostream &logfile)
{
	logfile << "i: " << i << ", j: " << j << endl;
	logfile.flush();
}

void Tracer::OutputProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	cout << (count*100)/(betaNumber+1) << '%'
		 << '\t' << timer.Elapsed();
}


void Tracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber, const string &path)
{
	string startTime = ctime(&m_startTime);
	string totalTime = timer.Elapsed();
	time_t end = timer.Stop();
	string endTime = ctime(&end);

	m_summary += "\nStart of calculation = " + startTime
			+ "End of calculation   = " + endTime
			+ "\nTotal time of calculation = " + totalTime
			+ "\nTotal number of body orientation = " + to_string(orNumber);

//	if (isNanOccured)
//	{
//		m_summary += "\n\nWARNING! NAN values occured. See 'log.txt'";
//	}

	ofstream out(path + "\\out.dat", ios::out);

	out << m_summary;
	out.close();

	cout << m_summary;
}

void Tracer::SetIsOutputGroups(bool value)
{
	isOutputGroups = value;
}

void Tracer::OutputStartTime(CalcTimer &timer)
{
	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;
}

void Tracer::SetHandler(Handler *handler)
{
	m_handler = handler;
	m_handler->SetScattering(m_scattering);
}

