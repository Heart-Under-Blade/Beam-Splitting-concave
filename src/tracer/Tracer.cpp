#include "Tracer.h"

#include <ostream>
#include <iostream>
#include <assert.h>
#include "common.h"
#include "macro.h"
//#ifdef _TRACK_ALLOW
//std::ofstream trackMapFile("tracks_deb.dat", std::ios::out);
//#endif

using namespace std;

LightTracer::LightTracer(Particle *particle, Scattering *scattering,
						 const string &resultFileName)
	: m_particle(particle),
	  m_resultDirName(resultFileName),
	  m_scattering(scattering)
{
}

LightTracer::~LightTracer()
{
}

void LightTracer::TraceFixed(const Angle3d &orientation)
{
	Angle3d orient = orientation.ToRadian();

	vector<Beam> outBeams;
	m_particle->Rotate(orient);
	m_scattering->ScatterLight(outBeams);
//	m_particle->Output();
	m_handler->HandleBeams(outBeams);
	outBeams.clear();

//	double D_tot = CalcTotalScatteringEnergy();

	m_handler->WriteMatricesToFile(m_resultDirName);
//	WriteStatisticsToFileGO(1, D_tot, 1, timer); // TODO: раскомментить
}

void LightTracer::TraceRandom(const AngleRange &/*betaRange*/,
							  const AngleRange &/*gammaRange*/)
{
}

void LightTracer::OutputOrientationToLog(int i, int j, ostream &logfile)
{
	logfile << "i: " << i << ", j: " << j << endl;
	logfile.flush();
}

void LightTracer::OutputProgress(int betaNumber, long long count, CalcTimer &timer)
{
	EraseConsoleLine(50);
	cout << (count*100)/(betaNumber+1) << '%'
		 << '\t' << timer.Elapsed();
}

void LightTracer::OutputStatisticsPO(CalcTimer &timer, long long orNumber, const string &path)
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

void LightTracer::SetIsOutputGroups(bool value)
{
	isOutputGroups = value;
}

void LightTracer::OutputStartTime(CalcTimer &timer)
{
	m_startTime = timer.Start();
	cout << "Started at " << ctime(&m_startTime) << endl;
}

void LightTracer::SetHandler(Handler *handler)
{
	m_handler = handler;
	m_handler->SetScattering(m_scattering);
}

