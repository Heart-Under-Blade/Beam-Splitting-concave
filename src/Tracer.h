#pragma once

#include "geometry_lib.h"
#include "PhysMtr.hpp"
#include "Tracing.h"
#include "CalcTimer.h"
#include "Mueller.hpp"
#include "BigInteger.hh"
#include "MullerMatrix.h"

class PointContribution;
class ScatteringFiles;

struct AngleRange
{
	double min;
	double max;
	int number;
	double norm;
	double step;

	AngleRange(double _min, double _max, int _number)
		: number(_number)
	{
		min = _min;
		max = _max;
		norm = max - min;
		step = norm/number;
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
	double m_outcomingEnergy;

	Light m_incidentLight;
	std::string m_resultDirName;
	double normIndex;
	double m_wavelength;
	Symmetry m_symmetry;
	std::string m_summary;
	time_t m_startTime;

	// REF: заменить
	bool isCalcOther = false;

	bool isNan = false;
	bool isOutputGroups = false;

private:
	std::vector<Arr2DC> J; // Jones matrices

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

	void AllocJ(std::vector<Arr2DC> &j, int m, int n, int size);
	void CleanJ(std::vector<Arr2DC> &j);
	void OutputToAllFile(std::ofstream &diffFile, std::ofstream &otherFile,
						 double degBeta, std::ofstream &allFile, Arr2D &all,
						 Arr2D &other);
	void SetIncidentLight(Particle *particle);
};
