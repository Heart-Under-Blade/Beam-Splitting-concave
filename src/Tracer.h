#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"
#include "CalcTimer.h"
#include "Mueller.hpp"

struct TrackGroup
{
	int groupID;
	long long int arr[1024];
	int size = 0;
	std::vector<std::vector<int>> tracks;
};

struct Tracks
{
	TrackGroup groups[32];
	int count = 0;

	int GetMaxGroupID() const // REF: вызывать в конструкторе
	{
		int maxGroupID = 0;

		for (int i = 0; i < count; ++i)
		{
			if (groups[i].groupID > maxGroupID)
			{
				maxGroupID = groups[i].groupID;
			}
		}

		return ++maxGroupID;
	}

	int GetGroupID(long long int trackID) const
	{
		for (int i = 0; i < count; ++i)
		{
			for (int j = 0; j < groups[i].size; ++j)
			{
				if (groups[i].arr[j] == trackID)
				{
					return groups[i].groupID;
				}
			}
		}

		if (count == 0)
		{
			return 0;
		}

		return -1;
	}
};

struct AngleRange
{
	int count;
	double norm;
	double step;

	AngleRange(double begin, double end, int count, double normCoef)
		: count(count)
	{
		norm = normCoef/count;
		step = DegToRad(end - begin)/count;
	}
};

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

class Tracer
{
public:
	Tracer(Tracing *tracing, const std::string resultFileName);

	void TraceIntervalGO(const AngleRange &betaR, const AngleRange &gammaR, int thetaNum);
	void TraceSingleOrGO(const double &beta, const double &gamma,
						 int thetaNum, const Tracks &tracks);

	void TraceIntervalPO(const AngleRange &betaR, const AngleRange &gammaR,
						 const Cone &bsCone, const Tracks &tracks, double wave);
	void TraceIntervalPO2(const AngleRange &betaR, const AngleRange &gammaR,
						 const Cone &bsCone, const Tracks &tracks, double wave);
	void TraceSingleOrPO(const double &beta, const double &gamma,
						 const Cone &bsCone, const Tracks &tracks, double wave);

private:
	Tracing *m_tracing;
	Arr2D m_muller;
	Point3f m_incidentDir;
	Point3f m_polarizationBasis;
	std::string m_resultFileName;
	double m_gammaNorm;
	std::vector<Arr2DC> J; // Jones matrices
	double sizeBin;

	matrix back;	///< Mueller matrix in backward direction
	matrix forw;	///< Mueller matrix in forward direction

	// light energy balance
	double m_incomingEnergy;
	double m_outcomingEnergy;
	time_t m_startTime;

private:
	void CleanJ(int maxGroupID, const Cone &bsCone);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta, const Tracks &tracks);
	void HandleBeamsPO(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, const Tracks &tracks);
	void HandleBeamsPO2(std::vector<Beam> &outBeams, const Cone &bsCone, double wavelength, int groupID);
	void SetJnRot(Beam &beam, const Point3f &T,
				  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot);
	void AddResultToSumMatrix(Arr2D &M_, int maxGroupID, const Cone &bsCone,
							  double norm = 1);
	void WriteSumMatrix(std::ofstream &outFile, const Arr2D &sum,
						const Cone &bsCone);
	void AddToSumMatrix(const Cone &bsCone, double norm, int q, Arr2D &M_);
	void PrintProgress(int betaNumber, long long count, CalcTimer &timer);
	void ExtractPeaksGO(int EDF, double NRM, int ThetaNumber);
	void WriteResultsToFileGO(int thetaNum, double NRM, const std::string &filename);
	void WriteStatisticsToFileGO(int orNumber, double D_tot, double NRM);
	std::string GetFileName(const std::string &filename);

	double CalcNorm(long long orNum);
	double CalcTotalScatteringEnergy(int thetaNum);
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AddToResultMullerGO(const Point3f &dir, matrix &bf, double area);
};
