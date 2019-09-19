#pragma once

#include "Beam.h"
#include "Scattering.h"
#include "PhysMtr.hpp"
#include "MullerMatrix.h"
#include "Tracks.h"
#include "ScatteringFiles.h"
#include "common.h"

class ScatteringSphere
{
public:
	ScatteringSphere(double radius, int phiCount, int thetaCount)
		: radius(radius), nAzimuth(phiCount), nZenith(thetaCount)
	{
		azinuthStep = M_2PI/phiCount;
		zenithStep = Orientation::DegToRad(radius)/thetaCount;
	}

	void ComputeSphereDirections(const Point3f &polarizationBasis)
	{
		double sinAz;
		double cosAz;

		double sinZen;
		double cosZen;

		for (int i = 0; i <= nAzimuth; ++i)
		{
			double p = i * azinuthStep;
			sincos(p, &sinAz, &cosAz);

			vf.push_back(Point3d(-sinAz ,cosAz ,0));
			directions.push_back(std::vector<Point3d>());

			for (int j = 0; j <= nZenith; ++j)
			{
				double t = j * zenithStep;
				sincos(t, &sinZen, &cosZen);

				Point3d dir(sinZen*cosAz, sinZen*sinAz, cosZen);
				directions[i].push_back(dir);
			}
		}

		vf.push_back(-polarizationBasis);
	}

public:
	double radius;
	int nAzimuth;
	int nZenith;
	double azinuthStep;
	double zenithStep;

	std::vector<Point3d> vf;
	std::vector<std::vector<Point3d>> directions;
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

class Handler
{
public:
	Handler(Scattering *scattering, double wavelength = 0);

	virtual void HandleBeams(std::vector<Beam> &beams);
	virtual void SetTracks(Tracks *tracks);
	Tracks *GetTracks() const;
	void SetScattering(Scattering *scattering);
	virtual void WriteMatricesToFile(std::string &destName);
	void EnableAbsorption(bool isNew);

	void SetNormIndex(double value);
	void SetSinZenith(double value);

	std::ofstream m_logFile;

	int m_nBadBeams;
	bool m_isBadBeam;

	void OutputPaths(BeamInfo &info, const OpticalPath &path);

	Beam m_startBeam;
	double m_sinZenith;
	double m_wavelength; // must be double type!!!
	Tracks *m_tracks;
	ScatteringSphere m_sphere; // back scattering conus

	virtual BeamInfo ComputeBeamInfo(Beam &beam);

protected:
	double BeamCrossSection(const Beam &beam) const;
	void ApplyAbsorption(BeamInfo &info);

	/**
	 * @brief Calculate the diffraction of beam in given direction
	 * @param beam beam of light
	 * @param direction direction for calculating diffraction
	 * @return Fresnel coefficient
	 */
	complex DiffractIncline(const BeamInfo &info, const Beam &beam,
							const Point3d &direction) const;

	complex DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
							   const Point3d &direction) const;

	Point3d ChangeCoordinateSystem(const Axes &axes, const Point3d& normal,
								   const Point3d& point) const;

	Point3d ChangeCoordinateSystem(const Point3d& normal,
								   const Point3d &point) const;

	void ComputeCoordinateSystemAxes(const Point3d& normal, Axes &axes) const;

	void ComputeLengthIndices(const Beam &beam, BeamInfo &info);

	void ComputeOpticalLengths(const Beam &beam, BeamInfo &info);

protected:
	Scattering *m_scattering;

	Particle *m_particle;
	bool m_hasAbsorption;
	double m_normIndex;
	double m_cAbs;
	complex m_cAbsExp;
	bool m_isNewAbs;
	std::ofstream m_absLogFile;

	complex m_ri;
	double m_riIm;

	double m_wavenumber;
	double m_wn2;

	complex m_complWave;
	complex m_invComplWave;

	double m_eps1;
	double m_eps2;
	double m_eps3;

	int count = 0;

private:
	void ExtropolateOpticalLenght(Beam &beam, const std::vector<int> &tr);
	complex ComputeAvgBeamEnergy(const Polygon &pol, const BeamInfo &info,
								 Point3d &p1, Point3d &p2,
								 double &p11, double &p12,
								 double &p21, double &p22,
								 const complex &c1, const complex &c2) const;
};
