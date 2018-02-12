#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"
#include "CalcTimer.h"
#include "Mueller.hpp"
#include "BigInteger.hh"
#include "MullerMatrix.h"

#define MAX_GROUP_NUM 1024

class PointContribution;
class ScatteringFiles;

struct ContributionGO
{
	ContributionGO()
		: scatMatrix(0, 0, 0, 0), back(4, 4), forw(4, 4)
	{
		scatMatrix = Arr2D(1, 180 + 1/*TODO: this is 'thetaNum',
				   should i do smth with it?*/, 4, 4);
		scatMatrix.ClearArr();
	}

	Arr2D scatMatrix;	///< Scattering matrix
	matrix back;		///< Mueller matrix in backward direction
	matrix forw;		///< Mueller matrix in forward direction
};

struct TrackGroup
{
	int groupID;
	BigInteger arr[MAX_GROUP_NUM];
	int size = 0;
	std::vector<std::vector<int>> tracks;

	std::string CreateGroupName() const
	{
		std::string subname;
		subname += "gr_" + std::to_string(groupID);
		return subname;
	}
};

class Tracks : public std::vector<TrackGroup>
{
public:
	int FindGroup(const BigInteger &trackID) const
	{
		for (size_t i = 0; i < size(); ++i)
		{
			for (int j = 0; j < (*this)[i].size; ++j)
			{
				if ((*this)[i].arr[j] == trackID)
				{
					return (*this)[i].groupID;
				}
			}
		}

		if (size() == 0)
		{
			return 0;
		}

		return -1;
	}

	static void RecoverTrack(const Beam &beam, int facetNum,
							 std::vector<int> &track)
	{
		int coef = facetNum + 1;
		std::vector<int> tmp_track;

		BigInteger tmpId = beam.id/coef;
		for (int i = 0; i <= beam.level; ++i)
		{
			int tmp = (tmpId%coef).toInt();
			tmpId -= tmp;
			tmpId /= coef;
			tmp -= 1;
			tmp_track.push_back(tmp);
		}

		for (int i = tmp_track.size()-1; i >= 0; --i)
		{
			track.push_back(tmp_track.at(i));
		}
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
	Tracer(Particle *particle, int reflNum, const std::string &resultFileName);
	~Tracer();

	void TraceRandomGO(int betaNumber, int gammaNumber, const Tracks &tracks);
	void TraceRandomGO(int betaNumber, int gammaNumber);
	void TraceFixedGO(const double &beta, const double &gamma,
					  const Tracks &tracks);

	void TraceRandomPO(int betaNumber, int gammaNumber, const Cone &bsCone,
					   const Tracks &tracks, double wave);

	void TraceRandomPO2(int betaNumber, int gammaNumber, const Cone &bsCone,
						const Tracks &tracks, double wave);
	void TraceFixedPO(const double &beta, const double &gamma,
					  const Cone &bsCone, const Tracks &tracks, double wave);

	void SetIsCalcOther(bool value); // REF: заменить
	void SetIsOutputGroups(bool value);// REF: заменить

	void OutputStatisticsPO(CalcTimer &timer, long long orNumber, const std::string &path);


protected:
	Tracing *m_tracing;
	double m_incomingEnergy;
	std::string m_resultDirName;
	double normIndex;
	double m_wavelength;
	Point3f m_polarizationBasis;
	Point3f m_incidentDir;

	// REF: заменить
	bool isCalcOther = false;

	bool isNan = false;
	bool isOutputGroups = false;

private:
	Symmetry m_symmetry;

	// result scattering martices
	ContributionGO m_totalMtrx;
	std::vector<ContributionGO> m_sepatateMatrices; // матрицы для вклада в отдельные траектории

	std::vector<Arr2DC> J; // Jones matrices

	// light energy balance
	double m_outcomingEnergy;
	time_t m_startTime;


	std::string m_statistics;

	bool isNanOccured = false;

protected:
	void OutputStartTime(CalcTimer &timer);
	void OutputProgress(int betaNumber, long long count, CalcTimer &timer);
	void OutputOrientationToLog(int i, int j, std::ostream &logfile);
	void CalcMultiplyOfJmatrix(const Beam &beam, const Point3f &T,
							   const Point3d &vf, const Point3d &vr,
							   double lng_proj0, matrixC &Jx);

private:
	void CleanJ(int size, const Cone &bsCone);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta);
	void HandleBeamsGO(std::vector<Beam> &outBeams, double beta, const Tracks &tracks);
	void HandleBeamsPO(std::vector<Beam> &outBeams, const Cone &bsCone,
					   const Tracks &tracks);
	void HandleBeamsPO2(std::vector<Beam> &outBeams, const Cone &bsCone, int groupID);

	void CalcJnRot(const Beam &beam, const Point3f &T,
				   const Point3d &vf, const Point3d &vr, matrixC &Jn_rot);
	void AddResultToMatrix(Arr2D &M, const Cone &bsCone, double norm = 1);
	void AddResultToMatrix(Arr2D &M, std::vector<Arr2DC> &j, double norm = 1);
	void AddResultToMatrices(std::vector<Arr2D> &M, const Cone &bsCone,
							 double norm = 1);
	void AddResultToMatrices(std::vector<Arr2D> &M);
	void WriteConusMatrices(std::ofstream &outFile, const Arr2D &sum,
						const Cone &bsCone);
	void AddToSumMatrix(const Cone &bsCone, double norm, int q, Arr2D &M_);
	void ExtractPeaksGO(int EDF, double NRM, ContributionGO &contr);
	void WriteResultsToFileGO(double NRM, const std::string &filename,
							  ContributionGO &contr);

	double CalcNorm(long long orNum);
	double CalcTotalScatteringEnergy();
	void RotateMuller(const Point3f &dir, matrix &bf);
	void AddToResultMullerGO(const Point3f &dir, matrix &bf, double area,
							 ContributionGO &contr);
	void WriteResultToSeparateFilesGO(double NRM, int EDF, const std::string &dir,
									  const Tracks &tracks);
	void AllocJ(std::vector<Arr2DC> &j, int m, int n, int size);
	void CleanJ(std::vector<Arr2DC> &j);
	void OutputStatisticsGO(int orNumber, double D_tot, double NRM,
						  CalcTimer &timer);
	void OutputToAllFile(std::ofstream &diffFile, std::ofstream &otherFile,
						 double degBeta, std::ofstream &allFile, Arr2D &all,
						 Arr2D &other);
};
