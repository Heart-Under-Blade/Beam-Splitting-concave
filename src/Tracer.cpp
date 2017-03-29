#include "Tracer.h"

#include <iostream>
#include "Mueller.hpp"
#include "CalcTimer.h"

using namespace std;

Tracer::Tracer(Tracing *tracing, const string resultFileName)
	: m_tracing(tracing),
	  m_mxd(0, 0, 0, 0),
	  m_incidentDir(0, 0, -1),
	  m_polarizationBasis(0, 1, 0),
	  m_resultFileName(resultFileName),
	  m_gammaNorm(3/M_PI)
{
}

void Tracer::TraceIntervalPO(const AngleInterval &betaI, const AngleInterval &gammaI,
							 const Cone &bsCone, const Tracks &tracks, double wave)
{
	CalcTimer timer;
	long long count = 0;

	Arr2D M(bsCone.phiCount+1, bsCone.thetaCount+1, 4, 4);
	ofstream outFile(m_resultFileName, ios::out);

	vector<Beam> outBeams;
	double beta, gamma;
	double step = betaI.GetStep();

	int maxGroupID = tracks.GetMaxGroupID();
	int halfGammaCount = gammaI.count/2;

	timer.Start();

	for (int i = 0; i <=betaI.count; ++i)
	{
		beta = step*i;

		for (int j = -halfGammaCount; j <= halfGammaCount; ++j)
		{
			gamma = j*gammaI.norm + M_PI/6;
//// DEB
//beta = DegToRad(32);
//gamma = DegToRad(30);
			m_tracing->SplitBeamByParticle(beta, gamma, outBeams);
			CleanJ(maxGroupID, bsCone);
			HandleBeamsPO(outBeams, bsCone, wave, tracks);
			outBeams.clear();
			AddResultToSumMatrix(M, maxGroupID, bsCone, gammaI);
		}

		WriteSumMatrix(outFile, M, bsCone);

		cout << '\r';

		for (int i = 0; i < 80; ++i)
		{
			cout << " ";
		}

		cout << '\r';

		cout << ((100*count)/betaI.count) << "% ";
	}

	outFile.close();
}

void Tracer::WriteSumMatrix(ofstream &outFile, const Arr2D &sum,
							const Cone &bsCone)
{
	outFile << to_string(bsCone.radius)	<< ' '
			<< to_string(bsCone.thetaCount) << ' '
			<< to_string(bsCone.phiCount+1);

	for (int t = 0; t <= bsCone.thetaCount; ++t)
	{
		double tt = (double)(t*bsCone.dTheta*180.0)/M_PI;

		for (int p = 0; p <= bsCone.phiCount; ++p)
		{
			double fi = -((double)p)*bsCone.dPhi;
			matrix m = sum(p ,t);
			double degPhi = RadToDeg(-fi);
			outFile << endl << tt << " " << degPhi << " ";
			outFile << m;
		}
	}
}

void Tracer::AddResultToSumMatrix(Arr2D &M_, int maxGroupID, const Cone &bsCone,
								  const AngleInterval &gammaI)
{
	double coef = gammaI.norm*m_gammaNorm;

	for (int q = 0; q < maxGroupID; ++q)
	{
		for (int t = 0; t <= bsCone.thetaCount; ++t)
		{
			for (int p = 0; p <= bsCone.phiCount; ++p)
			{
//				complex ee = J[q](p, t)[0][0]; // DEB
				matrix Mk = Mueller(J[q](p, t));
				M_.insert(p, t, coef*Mk);
			}
		}
	}
}

void Tracer::CleanJ(int maxGroupID, const Cone &bsCone)
{
	J.clear();
	Arr2DC tmp(bsCone.phiCount+1, bsCone.thetaCount+1, 2, 2);
	tmp.ClearArr();

	for(int q = 0; q < maxGroupID; q++)
	{
		J.push_back(tmp);
	}
}

void Tracer::TraceIntervalGO(const AngleInterval &betaI, const AngleInterval &gammaI)
{

}

void Tracer::TraceSingleOr(const double &beta, const double &gamma)
{
	vector<Beam> outBeams;
	m_tracing->SplitBeamByParticle(beta, gamma, outBeams);
}

void Tracer::SetJnRot(Beam &beam, const Point3f &T,
					  const Point3d &vf, const Point3d &vr, matrixC &Jn_rot)
{
	Point3f normal = beam.polygon.Normal();

	Point3d vt = CrossProductD(vf, vr);
	vt = vt/LengthD(vt);

	Point3f cpNT = CrossProduct(normal, T);
	Point3f cpNE = CrossProduct(normal, beam.e);

	Point3d cpNTd = Point3d(cpNT.cx, cpNT.cy, cpNT.cz);
	Point3d cpNEd = Point3d(cpNE.cx, cpNE.cy, cpNE.cz);

	Jn_rot[0][0] = -DotProductD(cpNTd, vf); // OPT: похоже на SetJMatrix
	Jn_rot[0][1] = -DotProductD(cpNEd, vf);
	Jn_rot[1][0] =  DotProductD(cpNTd, vt);
	Jn_rot[1][1] =  DotProductD(cpNEd, vt);
}

void Tracer::HandleBeamsPO(vector<Beam> &outBeams, const Cone &bsCone,
						   double wavelength, const Tracks &tracks)
{
	for (unsigned int i = 0; i < outBeams.size(); ++i)
	{
		Beam &beam = outBeams.at(i);

		double ctetta = DotProduct(beam.direction, -m_incidentDir);

		if (ctetta < 0.17364817766693034885171662676931)
		{	// отбрасываем пучки, которые далеко от конуса направления назад
//			continue;// DEB
		}

		int groupID = tracks.GetGroupID(beam.id);

		if (groupID < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentDir, m_polarizationBasis);

		Point3f center = beam.polygon.Center();
		double lng_proj0 = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f T = CrossProduct(beam.e, beam.direction);
		T = T/Length(T); // базис выходящего пучка

		for (int i = 0; i <= bsCone.phiCount; ++i)
		{
			for (int j = 0; j <= bsCone.thetaCount; ++j)
			{
				double f = i*bsCone.dPhi;
				double t = j*bsCone.dTheta;

				double sinT = sin(t), sinF = sin(f), cosF = cos(f);

				Point3d vr(sinT*cosF, sinT*sinF, cos(t));
				Point3d vf = (j == 0) ? -m_polarizationBasis
									  : Point3d(-sinF ,cosF ,0);
				matrixC Jn_rot(2, 2);
				SetJnRot(beam, T, vf, vr, Jn_rot);

				complex fn(0, 0);
				fn = beam.DiffractionIncline(vr, wavelength);

				double dp = DotProductD(vr, Point3d(center));
				complex tmp = exp_im(M_2PI*(lng_proj0-dp)/wavelength);
				matrixC fn_jn = beam.J * tmp;

				matrixC c = fn*Jn_rot*fn_jn;
				J[groupID].insert(i, j, c);
			}
		}
	}
}
