#pragma once

#include "Handler.h"

#define SPHERE_RING_NUM	180 // number of rings no the scattering sphere
#define BIN_SIZE M_PI/SPHERE_RING_NUM

class HandlerGO : public Handler
{
public:
	HandlerGO(Particle *particle, Light *incidentLight, double wavelength = 0);

	void SetTracks(Tracks *tracks) override;

	double ComputeTotalScatteringEnergy();
	void WriteLog(const std::string &str);

	void MultiplyMueller(const Beam &beam, matrix &m);

protected:
	ContributionGO m_totalContrib;	// result scattering martices contribution
	std::vector<ContributionGO> m_tracksContrib; // group contibution of beams

protected:
	matrix ComputeMueller(float zenAng, Beam &beam);
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AverageOverAlpha(int EDF, double norm, ContributionGO &contrib,
						  const std::string &destDir);

	void WriteToFile(ContributionGO &contrib, double norm,
					 const std::string &filename);

	double ComputeOpticalPathAbsorption(const Beam &beam);
	Point3f CalcK(std::vector<int> &tr);

private:
	void ExtractPeaks(double *b, double *f, double norm, const std::string &destDir);
};
