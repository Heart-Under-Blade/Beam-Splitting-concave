#include "HandlerPO.h"

#include "Mueller.hpp"

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, double wavelength)
	: Handler(particle, incidentLight, wavelength),
	  m_sphere(0.0, 0, 0)
{
	m_wavelength = 0.532;
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

	complex fresnel = (m_hasAbsorption)
			? DiffractInclineAbs(info, beam, direction)
			: DiffractIncline(info, beam, direction);
	matrixC difracted = fresnel*jones_rot*fnJones;
#ifdef _DEBUG // DEB
	complex ddd[4];
	ddd[0] = jones_rot[0][0];
	ddd[1] = jones_rot[0][1];
	ddd[2] = jones_rot[1][0];
	ddd[3] = jones_rot[1][1];

	complex qqq[4];
	qqq[0] = fnJones[0][0];
	qqq[1] = fnJones[0][1];
	qqq[2] = fnJones[1][0];
	qqq[3] = fnJones[1][1];

	complex bbb[4];
	bbb[0] = difracted[0][0];
	bbb[1] = difracted[0][1];
	bbb[2] = difracted[1][0];
	bbb[3] = difracted[1][1];
#endif
	return difracted;
}

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
								  const Vector3d &direction)
{
	double dp = DotProductD(direction, info.center);
	double arg = m_waveIndex*(info.projLenght - dp);
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
	if (m_hasAbsorption)
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

//	double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

//	if (path > DBL_EPSILON)
//	{
//		double abs = exp(m_cAbs*path);
//		beam.J *= abs;
//	}

	//	ExtropolateOpticalLenght(beam, tr);
}

void HandlerPO::HandleBeams(std::vector<Beam> &beams)
{
	CleanJ();
	int groupId = 0;

	for (Beam &beam : beams)
	{
#ifdef _DEBUG // DEB
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

		for (int i = 0; i <= m_sphere.nAzimuth; ++i)
		{
			for (int j = 0; j <= m_sphere.nZenith; ++j)
			{
#ifdef _DEBUG // DEB
				if (j == 160)
					int ff = 0;
#endif
				Point3d &dir = m_sphere.directions[i][j];
				Point3d &vf = (j == 0) ? m_sphere.vf.back() : m_sphere.vf[i];
				matrixC diffractedMatrix = ApplyDiffraction(beam, info, dir, vf);
				m_diffractedMatrices[groupId].insert(i, j, diffractedMatrix);
			}
		}
	}

	AddToMueller();
}
