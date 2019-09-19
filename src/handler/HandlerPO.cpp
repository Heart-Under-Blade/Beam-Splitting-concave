#include "HandlerPO.h"

#include <iostream>

#include "Mueller.hpp"\

HandlerPO::HandlerPO(Scattering *scattering, double wavelength)
	: Handler(scattering, wavelength)
{
	double path = m_scattering->GetFarFresnelZone() * 2;
	m_shadowBeamPhaseOffset = exp_im(m_wavenumber * path);
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

	if (m_tracks->empty())
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
		double tt = Orientation::RadToDeg(t*m_sphere.zenithStep);

		for (int p = 0; p <= m_sphere.nAzimuth; ++p)
		{
			double fi = -((double)p)*m_sphere.azinuthStep;
			double degPhi = Orientation::RadToDeg(-fi);
			outFile << std::endl << tt << " " << degPhi << " ";
			matrix m = M(p, t);
			outFile << m;
		}
	}
}

complex HandlerPO::ComputePhaseOffset(const BeamInfo &info,
									  const Vector3d &direction)
{
	double dp = Point3d::DotProduct(direction, info.center);
	double arg = m_wavenumber*(info.projLenght - dp);
	return exp_im(arg);
}

void HandlerPO::RotateJones(const Beam &beam, const BeamInfo &info,
							const Vector3d &vf, const Vector3d &direction,
							matrixC &matrix) const
{
	auto &dir = direction;

	Vector3d vt = Point3d::CrossProduct(vf, dir);
	vt = vt/Point3d::Length(vt);

	Vector3f NT = Point3f::CrossProduct(info.normal, info.beamBasis);
	Vector3f NE = Point3f::CrossProduct(info.normal, beam.polarizationBasis);

	Vector3d NTd = Vector3d(NT);
	Vector3d NPd = Vector3d(NE);

	Point3f DT = Point3f::CrossProduct(beam.direction, info.beamBasis);
	Point3f DP = Point3f::CrossProduct(beam.direction, beam.polarizationBasis);

	Point3d DTd = Point3d(DT);
	Point3d DPd = Point3d(DP);

	const Point3d nd = info.normal;

	Point3d cpT = Point3d::CrossProduct(dir, NTd)
			- Point3d::CrossProduct(dir, Point3d::CrossProduct(dir, Point3d::CrossProduct(nd, DTd)));
	Point3d cpP = Point3d::CrossProduct(dir, NPd)
			- Point3d::CrossProduct(dir, Point3d::CrossProduct(dir, Point3d::CrossProduct(nd, DPd)));

	matrix[0][0] = Point3d::DotProduct(cpT, vt)/2.0;
	matrix[0][1] = Point3d::DotProduct(cpP, vt)/2.0;
	matrix[1][0] = Point3d::DotProduct(cpT, vf)/2.0;
	matrix[1][1] = Point3d::DotProduct(cpP, vf)/2.0;
}

matrixC HandlerPO::ApplyDiffraction(const Beam &beam, const BeamInfo &info,
									const Vector3d &direction,
									const Vector3d &vf)
{
	complex phaseOffset = (info.isShadow) ? m_shadowBeamPhaseOffset
										  : ComputePhaseOffset(info, direction);
	matrixC jones_rot(2, 2);
	RotateJones(beam, info, vf, direction, jones_rot);

	complex fresnel = (m_hasAbsorption && !info.isShadow)
			? DiffractInclineAbs(info, beam, direction)
			: DiffractIncline(info, beam, direction);

	return fresnel*jones_rot*(beam.Jones * phaseOffset);
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
	m_sphere.ComputeSphereDirections(m_startBeam.polarizationBasis);
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams)
{
#ifdef _DEBUG // DEB
//	int cc = 0;
#endif
	CleanJ();
	int groupId = 0;

	for (Beam &beam : beams)
	{
		m_isBadBeam = false;
#ifdef _DEBUG // DEB
//		cc++;
//		if (cc == 330)
//			int ddddddd = 0;
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
					-m_startBeam.direction, m_startBeam.polarizationBasis);

		BeamInfo info = ComputeBeamInfo(beam);

		if (m_isBadBeam)
		{
			continue;
		}

		if (!m_isNewAbs && !info.isShadow && m_hasAbsorption && beam.actNo > 0)
		{
			ApplyAbsorption(info);
		}

		for (int i = 0; i <= m_sphere.nAzimuth; ++i)
		{
			for (int j = 0; j <= m_sphere.nZenith; ++j)
			{
				Point3d &dir = m_sphere.directions[i][j];
				Point3d &vf = (j == 0) ? m_sphere.vf.back() : m_sphere.vf[i];
				matrixC diffractedMatrix = ApplyDiffraction(beam, info, dir, vf);
#ifdef _DEBUG // DEB
				complex fff = diffractedMatrix[0][0];
				if (std::isnan(real(fff)))
					int fff = 0;
//					m_logFile << cc << std::endl;
#endif
				m_diffractedMatrices[groupId].insert(i, j, diffractedMatrix);
			}
		}
	}

	AddToMueller();
}
