#include "HandlerPO.h"

#include <iostream>

#include "Mueller.hpp"

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, double wavelength)
	: Handler(particle, incidentLight, wavelength)
{
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

void HandlerPO::WriteMatricesToFile(string &destName)
{
	ofstream outFile(destName, ios::app);

	outFile << std::to_string(m_sphere.radius) << ' '
			<< std::to_string(m_sphere.nZenith) << ' '
			<< std::to_string(m_sphere.nAzimuth+1);

	for (int t = 0; t <= m_sphere.nZenith; ++t)
	{
		double tt = Orientation::RadToDeg(t*m_sphere.zenithStep);

		for (int p = 0; p <= m_sphere.nAzimuth; ++p)
		{
			double fi = -((double)p)*m_sphere.azinuthStep;
<<<<<<< HEAD
			double degPhi = RadToDeg(-fi);
			outFile << std::endl << tt << " " << degPhi << " ";

=======
			double degPhi = Orientation::RadToDeg(-fi);
			outFile << endl << tt << " " << degPhi << " ";
>>>>>>> origin/refactor
			matrix m = M(p, t);
			outFile << m;
		}
	}
}

<<<<<<< HEAD
=======
matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
								  const Vector3d &direction)
{
	double dp = Point3d::DotProduct(direction, info.center);
	double arg = m_wavenumber*(info.projLenght - dp);
	return matrix*exp_im(arg);
}

>>>>>>> origin/refactor
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
	matrixC fnJones = ComputeFnJones(beam.Jones, info, direction);
	matrixC jones_rot(2, 2);
	RotateJones(beam, info, vf, direction, jones_rot);

<<<<<<< HEAD
	complex fresnel = (m_hasAbsorption && beam.lastFacetId != INT_MAX &&
			beam.nActs > 0)
=======
	complex fresnel = (m_hasAbsorption && !beam.IsShadow())
>>>>>>> origin/refactor
			? DiffractInclineAbs(info, beam, direction)
			: DiffractIncline(info, beam, direction);

	return fresnel*jones_rot*fnJones;
}

<<<<<<< HEAD
matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &matrix, const BeamInfo &info,
								  const Vector3d &direction)
{
	double dp = DotProductD(direction, info.center);
	double arg = m_wavenumber*(info.projLenght - dp);
	return matrix*exp_im(arg);
}

=======
>>>>>>> origin/refactor
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

<<<<<<< HEAD
=======
BeamInfo HandlerPO::ComputeBeamInfo(const Beam &beam)
{
	BeamInfo info;
	info.normal = beam.Normal();
	info.normald = Point3d(info.normal.coordinates[0],
			info.normal.coordinates[1], info.normal.coordinates[2]);

	info.order = Point3f::DotProduct(info.normal, beam.direction) > 0;

	if (!info.order)
	{
		info.normal = -info.normal;
		info.normald = -info.normald;
	}

	ComputeCoordinateSystemAxes(info.normald, info.horAxis, info.verAxis);

	info.center = beam.Center();
	info.projectedCenter = ChangeCoordinateSystem(info.horAxis, info.verAxis,
												  info.normald, info.center);
	if (m_hasAbsorption && !beam.IsShadow())
	{
		ComputeOpticalLengths(beam, info);
		ComputeLengthIndices(beam, info);
	}

	info.area = beam.Area();

	if (beam.IsShadow())
	{
		info.projLenght = 20000 + Point3d::DotProduct(info.center, beam.direction);
	}
	else
	{
		OpticalPath opticalPath = ComputeOpticalPath(beam);
		info.projLenght = opticalPath.GetTotal() + Point3d::DotProduct(info.center, beam.direction);
	}

	info.beamBasis = Point3f::CrossProduct(beam.polarizationBasis, beam.direction);
	info.beamBasis = info.beamBasis/Point3f::Length(info.beamBasis); // basis of beam

	return info;
}

>>>>>>> origin/refactor
void HandlerPO::SetScatteringSphere(const ScatteringSphere &grid)
{
	m_sphere = grid;
	M = Arr2D(m_sphere.nAzimuth+1, m_sphere.nZenith+1, 4, 4);

	m_sphere.ComputeSphereDirections(*m_incidentLight);
}

<<<<<<< HEAD
=======
void HandlerPO::ComputeOpticalLengths(const Beam &beam, BeamInfo &info)
{
	std::vector<int> tr;
	m_tracks->RecoverTrack(beam, tr);

	for (int i = 0; i < 3; ++i)
	{
		OpticalPath op = m_scattering->ComputeOpticalPath(beam, Point3f(info.center.x, info.center.y, info.center.z), tr);
		info.opticalLengths[i] = op.GetTotal();
	}
	//	ExtropolateOpticalLenght(beam, tr);
}

>>>>>>> origin/refactor
void HandlerPO::HandleBeams(std::vector<Beam> &beams)
{
#ifdef _DEBUG // DEB
//	int cc = 0;
#endif
	CleanJ();
	int groupId = 0;

	for (Beam &beam : beams)
	{
//		std::cout << cc++ << std::endl;
		m_isBadBeam = false;
#ifdef _DEBUG // DEB
//		cc++;
//		if (cc == 330)
//			int ddddddd = 0;
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

<<<<<<< HEAD
//		if (cc == 3)
//			std::cout << "-" << std::endl;

		if (m_isBadBeam)
		{
			continue;
		}

//		if (beam.lastFacetId != INT_MAX && beam.nActs > 0)
//		{
//			std::vector<int> tr;
//			Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

////			std::cout << tr.size() << std::endl;
//			double path = m_scattering->MeasureOpticalPath(beam, info.centerf, tr);

//			if (path > DBL_EPSILON)
//			{
//				double abs = exp(m_cAbs*path);
//				beam.J *= abs;
//			}
//		}
=======
		if (!beam.IsShadow())
		{
			OpticalPath path = ComputeOpticalPath(beam);

			if (path.GetTotal() > DBL_EPSILON)
			{
				double abs = exp(m_cAbs*path.GetTotal());
				beam.Jones *= abs;
			}
		}
>>>>>>> origin/refactor

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
					int fff = 0;
//					m_logFile << cc << std::endl;
#endif
				m_diffractedMatrices[groupId].insert(i, j, diffractedMatrix);
			}
		}
	}

	AddToMueller();
}
