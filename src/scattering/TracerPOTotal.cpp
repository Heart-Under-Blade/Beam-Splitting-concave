#include "TracerPOTotal.h"

#include <iostream>

using namespace std;

TracerPOTotal::TracerPOTotal(Particle *particle, Scattering *scattering,
							 const string &resultFileName)
	: TracerPO(particle, scattering, resultFileName)
{
}

void TracerPOTotal::TraceRandom(const OrientationRange &range)
{
	CalcTimer timer;
	long long count = 0;
	long long nOrientations = range.nZenith * range.nAzimuth;
	ofstream outFile(m_resultDirName + "_out.dat", ios::out);

	if (!outFile.is_open())
	{
		std::cerr << "Error! File \"" << m_resultDirName
				  << "\" was not opened. " << __FUNCTION__;
		throw std::exception();
	}

	vector<Beam> outBeams;
	Orientation angle;

	double normIndex = 2 * range.nAzimuth;
	m_handler->SetNormIndex(normIndex);

	++nOrientations;
	timer.Start();

	for (int i = 0; i <= range.nZenith; ++i)
	{
		angle.zenith = i*range.step.zenith;

		double sinZenith = (i == 0 || i == range.nZenith)
				? (1.0-cos(0.5*range.step.zenith))/normIndex
				: (cos((i-0.5)*range.step.zenith) -
				   cos((i+0.5)*range.step.zenith))/normIndex;

		m_handler->SetSinZenith(sinZenith);

		for (int j = 0; j < range.nAzimuth; ++j)
		{
			angle.azimuth = j*range.step.azimuth;

			m_particle->Rotate(angle);
			outBeams.push_back(m_scattering->SetShadowBeam());
			m_scattering->ScatterLight(outBeams);

			m_handler->HandleBeams(outBeams);
//			cout << "!" << endl;
			outBeams.clear();

			OutputProgress(nOrientations, count, timer);
			++count;
		}
	}

	std::string dir = CreateFolder(m_resultDirName);
	m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
	m_handler->WriteMatricesToFile(m_resultDirName);
	OutputStatisticsPO(timer, nOrientations, dir);
	outFile.close();
}
