#include "HandlerPO.h"

<<<<<<< HEAD
#include <iostream>

#include "Mueller.hpp"

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, double wavelength)
	: Handler(particle, incidentLight, wavelength)
{
//	m_wavelength = 0.532;
}

void HandlerPO::CleanJ()
{
	m_diffractedMatrices.clear();
	Arr2DC tmp(m_sphere.nAzimuth + 1, m_sphere.nZenith + 1, 2, 2);
	tmp.ClearArr();

	for (unsigned q = 0; q < m_tracks->size(); q++)
	{
		m_diffractedMatrices.push_back(tmp);
	}
}

void HandlerPO::WriteMatricesToFile(std::string &destName)
{
	std::ofstream outFile(destName, std::ios::app);

	outFile << std::to_string(m_sphere.radius) << ' '
			<< std::to_string(m_sphere.nZenith) << ' '
			<< std::to_string(m_sphere.nAzimuth+1);

	for (int t = 0; t <= m_sphere.nZenith; ++t)
	{
		double tt = RadToDeg(t*m_sphere.zenithStep);

		for (int p = 0; p <= m_sphere.nAzimuth; ++p)
		{
			double fi = -((double)p)*m_sphere.azinuthStep;
			double degPhi = RadToDeg(-fi);
			outFile << std::endl << tt << " " << degPhi << " ";
=======
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46

			matrix m = M(p, t);
			outFile << m;
		}
	}
}

<<<<<<< HEAD
void HandlerPO::RotateJones(const Beam &beam, const BeamInfo &info,
							const Vector3d &vf, const Vector3d &direction,
							matrixC &matrix) const
{
	auto &dir = direction;

	Vector3d vt = CrossProductD(vf, dir);
	vt = vt/LengthD(vt);

	Vector3f NT = CrossProduct(info.normal, info.beamBasis);
	Vector3f NE = CrossProduct(info.normal, beam.polarizationBasis);

	Vector3d NTd = Vector3d(NT.cx, NT.cy, NT.cz);
	Vector3d NPd = Vector3d(NE.cx, NE.cy, NE.cz);

	Point3f DT = CrossProduct(beam.direction, info.beamBasis);
	Point3f DP = CrossProduct(beam.direction, beam.polarizationBasis);

	Point3d DTd = Point3d(DT.cx, DT.cy, DT.cz);
	Point3d DPd = Point3d(DP.cx, DP.cy, DP.cz);

	const Point3d nd = info.normal;

	Point3d cpT = CrossProductD(dir, NTd)
			- CrossProductD(dir, CrossProductD(dir, CrossProductD(nd, DTd)));
	Point3d cpP = CrossProductD(dir, NPd)
			- CrossProductD(dir, CrossProductD(dir, CrossProductD(nd, DPd)));

	matrix[0][0] = DotProductD(cpT, vt)/2.0;
	matrix[0][1] = DotProductD(cpP, vt)/2.0;
	matrix[1][0] = DotProductD(cpT, vf)/2.0;
	matrix[1][1] = DotProductD(cpP, vf)/2.0;
}

matrixC HandlerPO::ApplyDiffraction(const Beam &beam, const BeamInfo &info,
									const Vector3d &direction,
									const Vector3d &vf)
{
	matrixC fnJones = ComputeFnJones(beam.J, info, direction);
	matrixC jones_rot(2, 2);
	RotateJones(beam, info, vf, direction, jones_rot);

	complex fresnel = (m_hasAbsorption && beam.lastFacetId != INT_MAX)
			? DiffractInclineAbs(info, beam, direction)
			: DiffractIncline(info, beam, direction);

	return fresnel*jones_rot*fnJones;
}

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
								  const Vector3d &direction)
{
	double dp = DotProductD(direction, info.center);
	double arg = m_wavenumber*(info.projLenght - dp);
	return matrix*exp_im(arg);
=======
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}

void HandlerPO::AddToMueller()
{
<<<<<<< HEAD
	for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
	{
		for (int t = 0; t <= m_sphere.nZenith; ++t)
		{
			for (int p = 0; p <= m_sphere.nAzimuth; ++p)
			{
				matrix m = Mueller(m_diffractedMatrices[q](p, t));
=======
	for (int q = 0; q < J.size(); ++q)
	{
		for (int t = 0; t <= m_conus.thetaCount; ++t)
		{
			for (int p = 0; p <= m_conus.phiCount; ++p)
			{
				matrix m = Mueller(J[q](p, t));
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
				m *= m_normIndex;
				M.insert(p, t, m);
			}
		}
	}
}

<<<<<<< HEAD
BeamInfo HandlerPO::ComputeBeamInfo(const Beam &beam)
{
	BeamInfo info;
	info.normal = beam.Normal();
	info.normald = Point3d(info.normal.cx, info.normal.cy, info.normal.cz);

	info.order = DotProduct(info.normal, beam.direction) > 0;

	if (!info.order)
	{
		info.normal = -info.normal;
		info.normald = -info.normald;
	}

	ComputeCoordinateSystemAxes(info.normald, info.horAxis, info.verAxis);

	info.center = beam.Center();
	info.projectedCenter = ChangeCoordinateSystem(info.horAxis, info.verAxis,
												  info.normald, info.center);
	if (m_hasAbsorption && beam.lastFacetId != INT_MAX)
	{
		ComputeOpticalLengths(beam, info);
		ComputeLengthIndices(beam, info);
	}

	info.area = beam.Area();

	info.projLenght = beam.opticalPath + DotProductD(info.center, beam.direction);

	info.beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
	info.beamBasis = info.beamBasis/Length(info.beamBasis); // basis of beam

	return info;
}

void HandlerPO::SetScatteringSphere(const ScatteringSphere &grid)
{
	m_sphere = grid;
	M = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);

	m_sphere.ComputeSphereDirections(*m_incidentLight);
}

void HandlerPO::ComputeOpticalLengths(const Beam &beam, BeamInfo &info)
{
	std::vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

	for (int i = 0; i < 3; ++i)
	{
		info.opticalLengths[i] = m_scattering->ComputeInternalOpticalPath(
					beam, beam.arr[i], tr);
	}
	//	ExtropolateOpticalLenght(beam, tr);
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams)
{
#ifdef _DEBUG // DEB
	int cc = 0;
#endif
	CleanJ();
	int groupId = 0;

	for (Beam &beam : beams)
	{
#ifdef _DEBUG // DEB
		cc++;
		if (cc == 330)
			int ddddddd = 0;
//		std::vector<int> tr;
//		Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
//		if (tr.size() == 2 && tr[0] == 2 && tr[1] == 4)
//			int fff = 0;
#endif
		if (m_tracks->shouldComputeTracksOnly)
		{
			groupId = m_tracks->FindGroupByTrackId(beam.id);

			if (groupId < 0)
			{
				continue;
			}
		}

		beam.polarizationBasis = beam.RotateSpherical(
					-m_incidentLight->direction,
					m_incidentLight->polarizationBasis);

		BeamInfo info = ComputeBeamInfo(beam);

		if (m_isBadBeam)
		{
			continue;
		}
//		std::cout << "2" << std::endl;
//		if (beam.lastFacetId != INT_MAX)
//		{
//			std::vector<int> tr;
//			Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

//			double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

//			if (path > DBL_EPSILON)
//			{
//				double abs = exp(m_cAbs*path);
//				beam.J *= abs;
//			}
//		}

		for (int i = 0; i <= m_sphere.nAzimuth; ++i)
		{
			for (int j = 0; j <= m_sphere.nZenith; ++j)
			{
				Point3d &dir = m_sphere.directions[i][j];
				Point3d &vf = (j == 0) ? m_sphere.vf.back() : m_sphere.vf[i];
				matrixC diffractedMatrix = ApplyDiffraction(beam, info, dir, vf);
#ifdef _DEBUG // DEB
				complex fff = diffractedMatrix[0][0];
				if (isnan(real(fff)))
					m_logFile << cc << std::endl;
#endif
				m_diffractedMatrices[groupId].insert(i, j, diffractedMatrix);
			}
		}
	}

	AddToMueller();
=======
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
>>>>>>> 03452a781c85ee0d91303dc91c948c61e251ec46
}
