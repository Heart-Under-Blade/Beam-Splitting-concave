#include "HandlerPO.h"

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

			matrix m = M(p, t);
			outFile << m;
		}
	}
}

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
}

void HandlerPO::AddToMueller()
{
	for (size_t q = 0; q < m_diffractedMatrices.size(); ++q)
	{
		for (int t = 0; t <= m_sphere.nZenith; ++t)
		{
			for (int p = 0; p <= m_sphere.nAzimuth; ++p)
			{
				matrix m = Mueller(m_diffractedMatrices[q](p, t));
				m *= m_normIndex;
				M.insert(p, t, m);
			}
		}
	}
}

void HandlerPO::SetScatteringSphere(const ScatteringSphere &grid)
{
	m_sphere = grid;
	M = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);

	m_sphere.ComputeSphereDirections(*m_incidentLight);
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
		m_isBadBeam = false;

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
}
