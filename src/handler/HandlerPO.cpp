#include "HandlerPO.h"

#include "Mueller.hpp"

using namespace std;

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, float wavelength)
	: Handler(particle, incidentLight, wavelength),
	  m_conus(0.0, 0, 0)
{
}

void HandlerPO::WriteMatricesToFile(string &destName)
{
	ofstream outFile(destName, ios::app);

	outFile << to_string(m_conus.radius) << ' '
			<< to_string(m_conus.thetaCount) << ' '
			<< to_string(m_conus.phiCount+1);

	for (int t = 0; t <= m_conus.thetaCount; ++t)
	{
		double tt = Orientation::RadToDeg(t*m_conus.dTheta);

		for (int p = 0; p <= m_conus.phiCount; ++p)
		{
			double fi = -((double)p)*m_conus.dPhi;
			double degPhi = Orientation::RadToDeg(-fi);
			outFile << endl << tt << " " << degPhi << " ";

			matrix m = M(p, t);
			outFile << m;
		}
	}
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams)
{
	CleanJ();

	for (Beam &beam : beams)
	{
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);

		Point3f center = beam.Center();

		auto path = ComputeOpticalPath(beam);
		double projLenght = path.GetTotal() + Point3f::DotProduct(center, beam.direction);

		Point3f beamBasis = Point3f::CrossProduct(beam.polarizationBasis, beam.direction);
		beamBasis = beamBasis/Point3f::Length(beamBasis); // basis of beam

		for (int i = 0; i <= m_conus.phiCount; ++i)
		{
			double p = i * m_conus.dPhi;
			double sinP = sin(p),
					cosP = cos(p);

			for (int j = 0; j <= m_conus.thetaCount; ++j)
			{	//
				double t = j * m_conus.dTheta;
				double sinT = sin(t);

				Point3d vr(sinT*cosP, sinT*sinP, cos(t));
				Point3d vf = (j == 0) ? -m_incidentLight->polarizationBasis
									  : Point3d(-sinP ,cosP ,0);
				// OPT: вышеописанные параметры можно вычислить один раз и занести в массив

				matrixC jones(0, 0);
				matrixC fnJones = ComputeFnJones(beam.Jones, center, vr, projLenght);
				ApplyDiffraction(beam, beamBasis, vf, vr, fnJones, jones);
				J[groupId].insert(i, j, jones);
			}
		}
	}

	AddToMueller();
}

void HandlerPO::SetScatteringConus(const Conus &conus)
{
	m_conus = conus;
	M = Arr2D(m_conus.phiCount + 1, m_conus.thetaCount + 1, 4, 4);
}

void HandlerPO::AddToMueller()
{
	for (int q = 0; q < J.size(); ++q)
	{
		for (int t = 0; t <= m_conus.thetaCount; ++t)
		{
			for (int p = 0; p <= m_conus.phiCount; ++p)
			{
				matrix m = Mueller(J[q](p, t));
				m *= m_normIndex;
				M.insert(p, t, m);
			}
		}
	}
}

void HandlerPO::setCon20(bool value)
{
	con20 = value;
}

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &jones, const Point3d &center,
								  const Vector3d &vr, double projLenght)
{
	double dp = Point3d::DotProduct(vr, center);
	double arg = M_2PI*(projLenght-dp)/m_wavelength;
	return jones * exp_im(arg);
}

void HandlerPO::ApplyDiffraction(const Beam &beam, const Point3f &beamBasis,
								 const Vector3d &vf, const Vector3d &vr,
								 const matrixC &fnJones, matrixC &jones)
{
	matrixC jones_rot(2, 2);
	RotateJones(beam, beamBasis, vf, vr, jones_rot);

	complex fresnel = beam.DiffractionIncline(vr, m_wavelength);

#ifdef _DEBUG // DEB
	if (isnan(real(fresnel)))
	{
		isNanOccured = isNan = true;
		return;
	}
#endif

	jones = fresnel*jones_rot*fnJones;
#ifdef _DEBUG // DEB
	Matrix2x2c mm(fnJones);
	Matrix2x2c jr(jones_rot);
	Matrix2x2c jo(jones);
	int ffff = 0;
#endif
}

void HandlerPO::RotateJones(const Beam &beam, const Vector3f &T,
							const Vector3d &vf, const Vector3d &vr, matrixC &J)
{
	Vector3f normal = beam.Normal();

	Vector3d vt = Point3d::CrossProduct(vf, vr);
	vt = vt/Point3d::Length(vt);

	Vector3f NT = Point3f::CrossProduct(normal, T);
	Vector3f NE = Point3f::CrossProduct(normal, beam.polarizationBasis);

	Vector3d NTd = Vector3d(NT.coordinates[0], NT.coordinates[1], NT.coordinates[2]);
	Vector3d NPd = Vector3d(NE.coordinates[0], NE.coordinates[1], NE.coordinates[2]);

//	J[0][0] = -DotProductD(NTd, vf);
//	J[0][1] = -DotProductD(NEd, vf);
//	J[1][0] =  DotProductD(NTd, vt);
//	J[1][1] =  DotProductD(NEd, vt);
	Point3f DT = Point3f::CrossProduct(beam.direction, T);
	Point3f DP = Point3f::CrossProduct(beam.direction, beam.polarizationBasis);

	Point3d DTd = Point3d(DT.coordinates[0], DT.coordinates[1], DT.coordinates[2]);
	Point3d DPd = Point3d(DP.coordinates[0], DP.coordinates[1], DP.coordinates[2]);

	Point3d nd = Point3d(normal.coordinates[0], normal.coordinates[1], normal.coordinates[2]);

	Point3d cpT = Point3d::CrossProduct(vr, NTd) - Point3d::CrossProduct(vr, Point3d::CrossProduct(vr, Point3d::CrossProduct(nd, DTd)));
	Point3d cpP = Point3d::CrossProduct(vr, NPd) - Point3d::CrossProduct(vr, Point3d::CrossProduct(vr, Point3d::CrossProduct(nd, DPd)));

	J[0][0] = Point3d::DotProduct(cpT, vt)/2.0;
	J[0][1] = Point3d::DotProduct(cpP, vt)/2.0;
	J[1][0] = Point3d::DotProduct(cpT, vf)/2.0;
	J[1][1] = Point3d::DotProduct(cpP, vf)/2.0;
}

void HandlerPO::CleanJ()
{
	J.clear();
	Arr2DC tmp(m_conus.phiCount + 1, m_conus.thetaCount + 1, 2, 2);
	tmp.ClearArr();

	for (int q = 0; q < m_tracks->size(); q++)
	{
		J.push_back(tmp);
	}
}
