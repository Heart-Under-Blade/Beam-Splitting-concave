#include "Handler.h"
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <string>

using namespace std;

Handler::Handler(Particle *particle, Light *incidentLight, float wavelength)
	: m_incidentLight(incidentLight),
	  m_particle(particle),
	  m_wavelength(wavelength),
	  m_hasAbsorbtion(false),
	  m_normIndex(1)
{
}

void Handler::HandleBeams(std::vector<Beam> &/*beams*/)
{
}

void Handler::SetTracks(Tracks *tracks)
{
	if (!m_tracks)
	{
		std::cerr << "Tracks are not set" << std::endl;
		throw std::exception();
	}

	m_tracks = tracks;
}

double Handler::BeamCrossSection(const Beam &beam)
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f normal = beam.facet->ex_normal; // normal of last facet of beam
	double cosA = Point3f::DotProduct(normal, beam.direction);
	double e = fabs(cosA);

	if (e < eps)
	{
		return 0;
	}

	double area = beam.Area();
	double len = Point3f::Length(normal);
	return (e*area)/len;
}

void Handler::WriteMatricesToFile(string &/*destName*/)
{
}

void Handler::SetNormIndex(double normIndex)
{
	m_normIndex = normIndex;
}

void Handler::SetAbsorbtionAccounting(bool value)
{
	m_hasAbsorbtion = value;
	m_cAbs = -M_2PI*imag(m_scattering->m_refractiveIndex)/m_wavelength;

	m_absLogFile.open("abslog1.txt", ios::out);
	m_absLogFile << setprecision(10);
	m_absLogFile << "No" << ' '
				 << "CtrPath" << ' '
				 << "AvgPath" << ' '
				 << "Max-Min" << ' '
				 << "Nact" << ' '
				 << "Npt" << ' '
				 << "Tr" << ' '
				 << endl;
}

void Handler::SetScattering(Scattering *scattering)
{
	m_scattering = scattering;
}

void Handler::OutputPaths(const Beam &beam, const OpticalPath &path)
{
	vector<int> track;
	m_tracks->RecoverTrack(beam, track);

	vector<double> ps;
	double sum = 0;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		OpticalPath p0 = m_scattering->ComputeOpticalPath(beam, beam.vertices[i],
														  track);
		ps.push_back(p0.internal);
		sum += p0.internal;
	}

	double maxPath = *std::max_element(ps.begin(), ps.end());
	double minPath = *std::min_element(ps.begin(), ps.end());
	double delta = fabs(path.internal - sum/ps.size());

	if (delta >= 10e-4)
	{
		m_absLogFile << ++count << ' '
					 << path.internal << ' '
					 << sum/ps.size() << ' '
					 << maxPath - minPath << ' '
					 << beam.actNo << ' '
					 << beam.nVertices << ' '
					 << Tracks::TrackToStr(track)
					 << endl;
	}
}

OpticalPath Handler::ComputeOpticalPath(const Beam &beam)
{
	vector<int> track;
	m_tracks->RecoverTrack(beam, track);
//	double opAbs = CalcOpticalPathAbsorption(beam);
	return m_scattering->ComputeOpticalPath(beam, beam.Center(), track);
}

void Handler::ApplyAbsorbtion(Beam &beam)
{
	auto path = ComputeOpticalPath(beam);

	if (path.internal > DBL_EPSILON)
	{
#ifdef _DEBUG // DEB
//		OutputPaths(beam, path);
//		if (fabs(path.GetTotal() - beam.opticalPath) >= 10e-4)
//			int ggg = 0;
#endif
		double abs = exp(m_cAbs*path.internal);
		beam.Jones *= abs;
	}
}
