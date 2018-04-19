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
struct Cone
{
	Cone(double radius, int phiCount, int thetaCount)
		: radius(radius), phiCount(phiCount), thetaCount(thetaCount)
	{
		dPhi = M_2PI/(phiCount+1);
		dTheta = DegToRad(radius)/thetaCount;
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
	PointContribution(size_t nGroups, double normIndex)
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

	void AddToGroup(const Matrix2x2c &jones, size_t groupId)
	{
		groupJones[groupId] += jones;
	}

	void SumGroupTotal()
	{
		for (size_t gr = 0; gr < m_nGroups; ++gr)
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

	const MuellerMatrix &GetGroupMueller(size_t groupID)
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

	size_t m_nGroups;
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


class Handler
{
public:
	Handler(Particle *particle, Light *incidentLight, float wavelength = 0);

	virtual void HandleBeams(std::vector<Beam> &beams);
	virtual void SetTracks(Tracks *tracks);
	void SetScattering(Scattering *scattering);
	virtual void WriteMatricesToFile(std::string &destName);

	void SetNormIndex(double normIndex);

	Light *m_incidentLight;

protected:
	void ApplyAbsorbtion(Beam &beam);

protected:
	Scattering *m_scattering;
	Particle *m_particle;
	Tracks *m_tracks;
	double m_cAbs;
	float m_wavelength;
	bool m_hasAbsorbtion;
	double m_normIndex;
	std::ofstream m_logFile;
};


class HandlerPO : public Handler
{
public:
	HandlerPO(Particle *particle, Light *incidentLight, float wavelength = 0);

	void HandleBeams(std::vector<Beam> &beams) override;
	void WriteMatricesToFile(std::string &destName) override;

	void SetScatteringConus(const Cone &conus);

protected:
	void MultiplyJones(const Beam &beam, const Point3f &T,
					   const Point3d &vf, const Point3d &vr,
					   double lng_proj0, matrixC &Jx);

	void RotateJones(const Beam &beam, const Point3f &T, const Point3d &vf,
					 const Point3d &vr, matrixC &J);

	void CleanJ();
	void AddToMueller();

protected:
	std::vector<Arr2DC> J;	// Jones matrices
	Arr2D M;				// Mueller matrices
	Cone m_conus;			// back scattering conus
	bool isNanOccured = false;
	bool isNan = false;
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
	void SetAbsorbtionAccounting(bool value);
	void WriteLog(const std::string &str);

	void MultiplyMueller(const Beam &beam, matrix &m);

protected:
	ContributionGO m_totalContrib;	// result scattering martices contribution
	std::vector<ContributionGO> m_tracksContrib; // group contibution of beams

	double m_sinAngle;

protected:
	double BeamCrossSection(const Beam &beam) const;
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
