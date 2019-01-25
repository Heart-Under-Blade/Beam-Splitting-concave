#pragma once

#include "Beam.h"
#include "Scattering.h"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
#include "Tracks.h"
#include "ScatteringFiles.h"

/**
 * @brief The Cone struct
 * Backscattering cone divided by cells
 */
struct Conus
{
	Conus(double radius, int phiCount, int thetaCount)
		: radius(radius), phiCount(phiCount), thetaCount(thetaCount)
	{
		dPhi = M_2PI/(phiCount+1);
		dTheta = Angle3d::DegToRad(radius)/thetaCount;
	}

	double radius;
	int phiCount;
	int thetaCount;
	double dPhi;
	double dTheta;
};


class PointContribution
{
public:
	PointContribution(int nGroups, double normIndex)
		: m_nGroups(nGroups),
		  m_normIndex(normIndex)
	{
		groupJones.resize(nGroups);
		groupMuellers.resize(nGroups);
		ResetJones();
	}

	void AddToMueller(const Matrix2x2c &jones)
	{
		MuellerMatrix m(jones);
		m *= m_normIndex;
		rest += m;
	}

	void AddToGroup(const Matrix2x2c &jones, int groupId)
	{
		groupJones[groupId] += jones;
	}

	void SumGroupTotal()
	{
		for (int gr = 0; gr < m_nGroups; ++gr)
		{
			MuellerMatrix m(groupJones[gr]);
			m *= m_normIndex;
			groupMuellers[gr] += m;
			groupTotal += m;
		}

		ResetJones();
	}

	void SumTotal()
	{
		total += groupTotal;
		total += rest;
	}

	const MuellerMatrix &GetGroupTotal() const
	{
		return groupTotal;
	}

	const MuellerMatrix &GetTotal() const
	{
		return total;
	}

	const MuellerMatrix &GetRest() const
	{
		return rest;
	}

	const MuellerMatrix &GetGroupMueller(int groupID)
	{
		return groupMuellers.at(groupID);
	}

	void Reset()
	{
		groupTotal.Reset();
		rest.Reset();
		total.Reset();

		for (MuellerMatrix &m : groupMuellers)
		{
			m.Reset();
		}
	}

private:
	std::vector<Matrix2x2c> groupJones;
	std::vector<MuellerMatrix> groupMuellers;

	MuellerMatrix groupTotal;
	MuellerMatrix rest;
	MuellerMatrix total;

	int m_nGroups;
	double m_normIndex;

	void ResetJones()
	{
		for (Matrix2x2c &j : groupJones)
		{
			j.Fill(0.f);
		}
	}
};


class ContributionGO
{
public:
	ContributionGO() : muellers(0, 0, 0, 0), back(4, 4), forward(4, 4)
	{
		muellers = Arr2D(1, 180 + 1/*TODO: this is 'thetaNum',
				   should i do smth with it?*/, 4, 4);
		muellers.ClearArr();
		back.Fill(0);
		forward.Fill(0);
	}

	void AddMueller(int angle, const matrix &m)
	{
		if (angle >= 180)
		{
			forward += m;
		}
		else if (angle <= 0)
		{
			back += m;
		}
		else
		{
			muellers.insert(0, angle, m);
		}
	}

	Arr2D muellers;		///< Scattering matrices
	matrix back;		///< Mueller matrix in backward direction
	matrix forward;		///< Mueller matrix in forward direction
};

class Handler // REF: rename to "Detector"
{
public:
	Handler(Particle *particle, Light *incidentLight, float wavelength = 0);

	virtual void HandleBeams(std::vector<Beam> &beams);
	virtual void SetTracks(Tracks *tracks);
	void SetScattering(Scattering *scattering);
	virtual void WriteMatricesToFile(std::string &destName);
	void SetAbsorbtionAccounting(bool value);
	void SetNormIndex(double normIndex);

	Light *m_incidentLight;
	Tracks *m_tracks;
	float m_wavelength;

	static double BeamCrossSection(const Beam &beam);

protected:
	void ApplyAbsorbtion(Beam &beam);
	void OutputPaths(const Beam &beam, const OpticalPath &path);
	OpticalPath ComputeOpticalPath(const Beam &beam);

protected:
	Scattering *m_scattering;

	Particle *m_particle;
	bool m_hasAbsorbtion;
	double m_normIndex;
	std::ofstream m_logFile;
	std::ofstream m_absLogFile;
	double m_cAbs;
	long long count = 0;
};


class HandlerPO : public Handler
{
public:
	HandlerPO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;

	void SetScatteringConus(const Conus &conus);

	void setCon20(bool value);

protected:
	void ApplyDiffraction(const Beam &beam, const Point3f &beamBasis,
						  const Vector3d &vf, const Vector3d &vr,
						  const matrixC &fnJones, matrixC &jones);

	void RotateJones(const Beam &beam, const Vector3f &T,
					 const Vector3d &vf, const Vector3d &vr, matrixC &J);

	void CleanJ();
	void AddToMueller();
	matrixC ComputeFnJones(const Matrix2x2c &jones, const Point3d &center,
						   const Vector3d &vr, double projLenght);

protected:
	std::vector<Arr2DC> J;	// Jones matrices
	Arr2D M;				// Mueller matrices
	Conus m_conus;			// back scattering conus
	bool isNanOccured = false;
	bool isNan = false;
	bool con20 = true;
};

class HandlerBackScatterPoint : public HandlerPO
{
public:
	HandlerBackScatterPoint(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void SetTracks(Tracks *tracks) override;

	void OutputContribution(ScatteringFiles &files, double angle, double energy,
							bool isOutputGroups, std::string prefix = "");

private:
	PointContribution *originContrib;
	PointContribution *correctedContrib;
};


class HandlerGO : public Handler
{
public:
	HandlerGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void SetTracks(Tracks *tracks) override;

	double ComputeTotalScatteringEnergy();
	void WriteLog(const std::string &str);

	void MultiplyMueller(const Beam &beam, matrix &m);

protected:
	ContributionGO m_totalContrib;	// result scattering martices contribution
	std::vector<ContributionGO> m_tracksContrib; // group contibution of beams

	double m_sinAngle;

protected:
	matrix ComputeMueller(int zenAng, Beam &beam);
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AverageOverAlpha(int EDF, double norm, ContributionGO &contrib);

	void WriteToFile(ContributionGO &contrib, double norm,
					 const std::string &filename);

	double ComputeOpticalPathAbsorption(const Beam &beam);
	Point3f CalcK(std::vector<int> &tr);

private:
	void ExtractPeaks(double *b, double *f, double norm);
};


class HandlerTotalGO : public HandlerGO
{
public:
	HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;

private:
	double ComputeInternalOpticalPath(const Beam &b, const std::vector<int> &tr);
};


class HandlerTracksGO : public HandlerGO
{
public:
	HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;
};
