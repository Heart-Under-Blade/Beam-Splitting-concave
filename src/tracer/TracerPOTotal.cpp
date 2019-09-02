#include "TracerPOTotal.h"

#include "handler/HandlerPO.h"
#include "Mueller.hpp"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
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
<<<<<<< HEAD:src/tracer/TracerPOTotal.cpp
	long long nOrientations = (betaRange.number + 1) * gammaRange.number;
=======
	long long nOrientations = range.nZenith * range.nAzimuth;
>>>>>>> origin/refactor:src/scattering/TracerPOTotal.cpp
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

<<<<<<< HEAD:src/tracer/TracerPOTotal.cpp
	for (int i = 0; i <= betaRange.number; ++i)
//	for (int i = 6; i <= 6/*betaRange.number*/; ++i)
=======
	for (int i = 0; i <= range.nZenith; ++i)
>>>>>>> origin/refactor:src/scattering/TracerPOTotal.cpp
	{
		angle.zenith = i*range.step.zenith;

<<<<<<< HEAD:src/tracer/TracerPOTotal.cpp
//#ifdef _DEBUG // DEB
//		beta = DegToRad(170);
//#endif
		double sinZenith = (i == 0 || i == betaRange.number)
				? (1.0-cos(0.5*betaRange.step))/normIndex
				: (cos((i-0.5)*betaRange.step) -
				   cos((i+0.5)*betaRange.step))/normIndex;
=======
		double sinZenith = (i == 0 || i == range.nZenith)
				? (1.0-cos(0.5*range.step.zenith))/normIndex
				: (cos((i-0.5)*range.step.zenith) -
				   cos((i+0.5)*range.step.zenith))/normIndex;
>>>>>>> origin/refactor:src/scattering/TracerPOTotal.cpp

//		std::cout << sinZenith << std::endl;
		m_handler->SetSinZenith(sinZenith);

<<<<<<< HEAD:src/tracer/TracerPOTotal.cpp
		for (int j = 0; j < gammaRange.number; ++j)
//		for (int j = 65; j < 67/*gammaRange.number*/; ++j)
		{
			gamma = j*gammaRange.step;
//			m_particle->Rotate(179.34, 37, 0);
			m_particle->Rotate(M_PI-beta, M_PI+gamma, 0);
			m_scattering->ExtractShadowBeam(outBeams);
			m_scattering->ScatterLight(outBeams);
//			std::cout << "0" << std::endl;
=======
		for (int j = 0; j < range.nAzimuth; ++j)
		{
			angle.azimuth = j*range.step.azimuth;

			m_particle->Rotate(angle);
			outBeams.push_back(m_scattering->SetShadowBeam());
			m_scattering->ScatterLight(outBeams);
>>>>>>> origin/refactor:src/scattering/TracerPOTotal.cpp

			m_handler->HandleBeams(outBeams);
//			std::cout << "exit" << std::endl;

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
