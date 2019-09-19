#include "TracerPOTotal.h"

#include <iostream>

#include "HandlerPO.h"
#include "Mueller.hpp"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"

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

			m_particle->Rotate(Orientation(M_PI-angle.zenith, M_PI+angle.azimuth));
			outBeams.push_back(m_scattering->ExtractShadowBeam());
			m_scattering->Scatter(&outBeams);

			m_handler->HandleBeams(outBeams);

			outBeams.clear();

			OutputProgress(nOrientations, count, timer);
			++count;

#ifdef _DEBUG // DEB
//			if (i == 9)
//			{
				double m = static_cast<HandlerPO*>(m_handler)->M(0,100)[0][0];
				m_handler->m_logFile << i << ' ' << j << ' ' << m << endl;
//			}
#endif
		}
	}

	std::string dir = CreateFolder(m_resultDirName);
	m_resultDirName = dir + m_resultDirName + '\\' + m_resultDirName;
	m_handler->WriteMatricesToFile(m_resultDirName);

//	std::ofstream outFile1(m_resultDirName + ".dat", std::ios::out);
//	for (int i = 0; i <= m_handler->m_sphere.nAzimuth; ++i)
//	{
//		for (int j = 0; j <= m_handler->m_sphere.nZenith; ++j)
//		{
//			matrix m(4,4);
//			m.Fill(0);
//			m = static_cast<HandlerPO*>(m_handler)->M(i,j);
//			outFile1 << i << " " << j << " ";
//			outFile1 << m << endl;
//		}
//	}

	OutputLogPO(timer, nOrientations, dir);
	outFile.close();
}
