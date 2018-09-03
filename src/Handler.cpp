#include "Handler.h"
#include "Mueller.hpp"
#include <iostream>
#include <limits>
#include <iomanip>

#define SPHERE_RING_NUM		180		// number of rings no the scattering sphere
#define BIN_SIZE			M_PI/SPHERE_RING_NUM
#define BEAM_DIR_LIM		0.9396

using namespace std;

Handler::Handler(Particle *particle, Light *incidentLight, float wavelength)
	: m_incidentLight(incidentLight),
	  m_particle(particle),
	  m_wavelength(wavelength),
	  m_hasAbsorbtion(false),
	  m_normIndex(1)
{
}

void Handler::HandleBeams(std::vector<Beam> &/*beams*/)
{
}

void Handler::SetTracks(Tracks *tracks)
{
	if (!m_tracks)
	{
		std::cerr << "Tracks are not set" << std::endl;
		throw std::exception();
	}

	m_tracks = tracks;
}

HandlerGO::HandlerGO(Particle *particle, Light *incidentLight, float wavelength)
	: Handler(particle, incidentLight, wavelength)
{
	m_logFile.open("log1.txt", ios::out);
	m_logFile << setprecision(numeric_limits<long double>::digits10 + 1);
}

void HandlerGO::ExtractPeaks(double *b, double *f, double norm)
{
	std::ofstream bck("back.dat", std::ios::out);
	std::ofstream frw("forward.dat", std::ios::out);
	frw << "M11 M22/M11 M33/M11 M44/M11";
	bck << "M11 M22/M11 M33/M11 M44/M11";

	if (f[0] <= DBL_EPSILON)
	{
		frw << "\n0 0 0 0";
	}
	else
	{
		frw << "\n" << f[0]*norm
			<< " " << f[1]/f[0]
			<< " " << f[1]/f[0]
			<< " " << f[2]/f[0];
	}

	if (b[0] <= DBL_EPSILON)
	{
		bck << "\n0 0 0 0";
	}
	else
	{
		bck << "\n" << b[0]*norm
			<< " " << b[1]/b[0]
			<< " " << -b[1]/b[0]
			<< " " << b[2]/b[0];
	}

	bck.close();
	frw.close();
}

void HandlerGO::AverageOverAlpha(int EDF, double norm, ContributionGO &contrib)
{
	//Analytical averaging over alpha angle
	double b[3], f[3];
	b[0] =  contrib.back[0][0];
	b[1] = (contrib.back[1][1] - contrib.back[2][2])/2.0;
	b[2] =  contrib.back[3][3];

	f[0] =  contrib.forward[0][0];
	f[1] = (contrib.forward[1][1] + contrib.forward[2][2])/2.0;
	f[2] =  contrib.forward[3][3];

	// Extracting the forward and backward peak in a separate file if needed
	if (EDF)
	{
		ExtractPeaks(b, f, norm);
	}
	else
	{
		contrib.muellers(0,SPHERE_RING_NUM,0,0) += f[0];
		contrib.muellers(0,0,0,0) += b[0];
		contrib.muellers(0,SPHERE_RING_NUM,1,1) += f[1];
		contrib.muellers(0,0,1,1) += b[1];
		contrib.muellers(0,SPHERE_RING_NUM,2,2) += f[1];
		contrib.muellers(0,0,2,2) -= b[1];
		contrib.muellers(0,SPHERE_RING_NUM,3,3) += f[2];
		contrib.muellers(0,0,3,3) += b[2];
	}
}

double HandlerGO::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f normal = m_particle->facets[beam.lastFacetId].ex_normal; // normal of last facet of beam
	double cosA = DotProduct(normal, beam.direction);
	double e = fabs(cosA);

	if (e < eps)
	{
		return 0;
	}

	double area = beam.Area();
	double len = Length(normal);
	return (e*area)/len;
}

void HandlerGO::MultiplyMueller(const Beam &beam, matrix &m)
{
	double cross = BeamCrossSection(beam);
	double area = cross*m_sinAngle;
	m *= area;
}

matrix HandlerGO::ComputeMueller(int zenAng, Beam &beam)
{
	matrix m = Mueller(beam.J);

	if (zenAng < 180 && zenAng > 0)
	{
		const float &y = beam.direction.cy;

		if (y*y > DBL_EPSILON)
		{	// rotate the Mueller matrix of the beam to appropriate coordinate system
			RotateMuller(beam.direction, m);
		}

#ifdef _CALC_AREA_CONTIBUTION_ONLY
		bf = matrix(4,4);
		bf.Identity();
#endif
	}

	MultiplyMueller(beam, m);
	return m;
}

void HandlerGO::RotateMuller(const Point3f &dir, matrix &bf)
{
	const float &x = dir.cx;
	const float &y = dir.cy;

	double tmp = y*y;
	tmp = acos(x/sqrt(x*x+tmp));

	if (y < 0)
	{
		tmp = M_2PI-tmp;
	}

	tmp *= -2.0;
	RightRotateMueller(bf, cos(tmp), sin(tmp));
}

void HandlerGO::WriteToFile(ContributionGO &contrib, double norm,
							const std::string &filename)
{
	string name = CreateUniqueFileName(filename);
	ofstream allFile(name, std::ios::out);

	allFile << "tetta M11 M12/M11 M21/M11 M22/M11 M33/M11 M34/M11 M43/M11 M44/M11";

	for (int j = SPHERE_RING_NUM; j >= 0; j--)
	{
		double tmp0 = 180.0/SPHERE_RING_NUM*(SPHERE_RING_NUM-j);
		double tmp1 = (j == 0) ? -(0.25*180.0)/SPHERE_RING_NUM : 0;
		double tmp2 = (j == (int)SPHERE_RING_NUM) ? (0.25*180.0)/SPHERE_RING_NUM : 0;

		// Special case in first and last step
		allFile << '\n' << tmp0 + tmp1 + tmp2;

		double sn = (j == 0 || j == (int)SPHERE_RING_NUM)
				? 1-cos(BIN_SIZE/2.0)
				: (cos((j-0.5)*BIN_SIZE)-cos((j+0.5)*BIN_SIZE));

		matrix bf = contrib.muellers(0, j);

		if (bf[0][0] <= DBL_EPSILON)
		{
			allFile << " 0 0 0 0 0 0 0 0";
		}
		else
		{
			allFile << ' ' << bf[0][0]*norm/(2.0*M_PI*sn)
					<< ' ' << bf[0][1]/bf[0][0]
					<< ' ' << bf[1][0]/bf[0][0]
					<< ' ' << bf[1][1]/bf[0][0]
					<< ' ' << bf[2][2]/bf[0][0]
					<< ' ' << bf[2][3]/bf[0][0]
					<< ' ' << bf[3][2]/bf[0][0]
					<< ' ' << bf[3][3]/bf[0][0];
		}
	}

	allFile.close();
}

Point3f HandlerGO::CalcK(vector<int> &tr)
{	// OPT: сделать из переменных ссылки
	Point3f k, tmp;
	Point3f n1 = m_particle->facets[tr[0]].in_normal;
	Point3f nq = m_particle->facets[tr[tr.size()-1]].in_normal;
	CrossProduct(nq, n1, tmp);
	CrossProduct(tmp, nq, k);

	for (int i = tr.size()-2; i > 0; --i)
	{
		Point3f ni = m_particle->facets[tr[i]].in_normal;
		k = k - ni*2*fabs(DotProduct(ni, k));
	}

	Normalize(k);
	return k;
}

double HandlerGO::ComputeOpticalPathAbsorption(const Beam &beam)
{	// OPT: вынести переменные из цикла
	double opticalPath = 0;

	vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

	Point3f k = CalcK(tr);

	Point3f n1 = m_particle->facets[tr[0]].in_normal;

	for (int i = 0; i < beam.size; ++i)
	{
		double delta = Length(beam.Center() - beam.arr[i])/Length(k);
		opticalPath += (delta*DotProduct(k, n1))/DotProduct(beam.direction, n1);
	}

	opticalPath /= beam.size;

	return opticalPath;
}

void Handler::WriteMatricesToFile(string &/*destName*/)
{
}

void Handler::SetNormIndex(double normIndex)
{
	m_normIndex = normIndex;
}

double HandlerGO::ComputeTotalScatteringEnergy()
{
	double D_tot = m_totalContrib.back[0][0] + m_totalContrib.forward[0][0];

	for (int i = 0; i <= SPHERE_RING_NUM; ++i)
	{
		D_tot += m_totalContrib.muellers(0, i, 0, 0);
	}

	return D_tot * m_normIndex;
}

void Handler::SetAbsorbtionAccounting(bool value)
{
	m_hasAbsorbtion = value;
	m_cAbs = -M_2PI*imag(m_particle->GetRefractiveIndex())/m_wavelength;
}

void HandlerGO::WriteLog(const string &str)
{
	m_logFile << str;
}

void Handler::SetScattering(Scattering *scattering)
{
	m_scattering = scattering;
}

HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void Handler::ApplyAbsorbtion(Beam &beam)
{
	vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

//	double opAbs = CalcOpticalPathAbsorption(beam);
#ifdef _DEBUG // DEB
	if (beam.id == 407490777/*415047647*/)
		int fff = 0;
#endif
	double path = m_scattering->ComputeInternalOpticalPath(beam, tr);

#ifdef _DEBUG // DEB
	double ddd = fabs(path - beam.opticalPath);
<<<<<<< HEAD
	m_logFile << ddd << " ";
	for (int t : tr) m_logFile << t << " ";
	m_logFile << endl;
=======
	m_logFile << ddd << endl;
>>>>>>> develop
	if (fabs(path - beam.opticalPath) >= 10e-4)
		int ggg = 0;
#endif
	if (path > DBL_EPSILON)
	{
		double abs = exp(m_cAbs*path);
		beam.J *= abs;
	}
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinAngle = sin(m_particle->rotAngle.beta);

	for (int i = 0; i < beams.size(); ++i)
	{
		Beam &beam = beams[i];
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);
		// absorbtion
<<<<<<< HEAD
		if (m_hasAbsorbtion && beam.act > 0)
=======
		if (m_hasAbsorbtion && beam.nActs > 0)
>>>>>>> develop
		{
			ApplyAbsorbtion(beam);
		}

		const float &z = beam.direction.cz;
		int zenith = round((acos(z)*SPHERE_RING_NUM)/M_PI);
		matrix m = ComputeMueller(zenith, beam);

		m_totalContrib.AddMueller(zenith, m);
	}
#ifdef _DEBUG // DEB
	double dd = m_totalContrib.muellers(0,0,0,0);
	m_logFile << dd;
	int ddd = 0;
#endif
}

HandlerTracksGO::HandlerTracksGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTracksGO::HandleBeams(std::vector<Beam> &beams)
{
	m_sinAngle = sin(m_particle->rotAngle.beta);

	for (Beam &beam : beams)
	{
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId >= 0)
		{
			beam.RotateSpherical(-m_incidentLight->direction,
								 m_incidentLight->polarizationBasis);

			const float &z = beam.direction.cz;
			int zenith = round((acos(z)*SPHERE_RING_NUM)/M_PI);
			matrix m = ComputeMueller(zenith, beam);

			m_totalContrib.AddMueller(zenith, m);
			m_tracksContrib[groupId].AddMueller(zenith, m);
		}
	}
}

void HandlerTracksGO::WriteMatricesToFile(string &destName)
{
	string dir = CreateFolder(destName);

	for (size_t i = 0; i < m_tracksContrib.size(); ++i)
	{
		if ((*m_tracks)[i].size != 0)
		{
			string subname = (*m_tracks)[i].CreateGroupName();
			AverageOverAlpha(0, m_normIndex, m_tracksContrib[i]);
			WriteToFile(m_tracksContrib[i], m_normIndex, dir + subname);
		}
	}
}

void HandlerTotalGO::WriteMatricesToFile(string &destName)
{
	destName += "_all";
	AverageOverAlpha(0, m_normIndex, m_totalContrib);
	WriteToFile(m_totalContrib, m_normIndex, destName);
}

HandlerPO::HandlerPO(Particle *particle, Light *incidentLight, float wavelength)
	: Handler(particle, incidentLight, wavelength),
	  m_conus(0.0, 0, 0)
{
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
		double projLenght = beam.opticalPath + DotProduct(center, beam.direction);

		Point3f beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
		beamBasis = beamBasis/Length(beamBasis); // basis of beam

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
				matrixC fnJones = ComputeFnJones(beam.J, center, vr, projLenght);
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
	for (size_t q = 0; q < J.size(); ++q)
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

matrixC HandlerPO::ComputeFnJones(const Matrix2x2c &jones, const Point3d &center,
								  const Vector3d &vr, double projLenght)
{
	double dp = DotProductD(vr, center);
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

	Vector3d vt = CrossProductD(vf, vr);
	vt = vt/LengthD(vt);

	Vector3f NT = CrossProduct(normal, T);
	Vector3f NE = CrossProduct(normal, beam.polarizationBasis);

	Vector3d NTd = Vector3d(NT.cx, NT.cy, NT.cz);
	Vector3d NPd = Vector3d(NE.cx, NE.cy, NE.cz);

//	J[0][0] = -DotProductD(NTd, vf);
//	J[0][1] = -DotProductD(NEd, vf);
//	J[1][0] =  DotProductD(NTd, vt);
//	J[1][1] =  DotProductD(NEd, vt);
	Point3f DT = CrossProduct(beam.direction, T);
	Point3f DP = CrossProduct(beam.direction, beam.polarizationBasis);

	Point3d DTd = Point3d(DT.cx, DT.cy, DT.cz);
	Point3d DPd = Point3d(DP.cx, DP.cy, DP.cz);

	Point3d nd = Point3d(normal.cx, normal.cy, normal.cz);

	Point3d cpT = CrossProductD(vr, NTd) - CrossProductD(vr, CrossProductD(vr, CrossProductD(nd, DTd)));
	Point3d cpP = CrossProductD(vr, NPd) - CrossProductD(vr, CrossProductD(vr, CrossProductD(nd, DPd)));

	J[0][0] = DotProductD(cpT, vt)/2.0;
	J[0][1] = DotProductD(cpP, vt)/2.0;
	J[1][0] = DotProductD(cpT, vf)/2.0;
	J[1][1] = DotProductD(cpP, vf)/2.0;
}

void HandlerPO::CleanJ()
{
	J.clear();
	Arr2DC tmp(m_conus.phiCount + 1, m_conus.thetaCount + 1, 2, 2);
	tmp.ClearArr();

	for (unsigned q = 0; q < m_tracks->size(); q++)
	{
		J.push_back(tmp);
	}
}

void HandlerGO::SetTracks(Tracks *tracks)
{
	Handler::SetTracks(tracks);
	m_tracksContrib.resize(m_tracks->size());
}

void HandlerPO::WriteMatricesToFile(string &destName)
{
	ofstream outFile(destName, ios::app);

	outFile << to_string(m_conus.radius) << ' '
			<< to_string(m_conus.thetaCount) << ' '
			<< to_string(m_conus.phiCount+1);

	for (int t = 0; t <= m_conus.thetaCount; ++t)
	{
		double tt = RadToDeg(t*m_conus.dTheta);

		for (int p = 0; p <= m_conus.phiCount; ++p)
		{
			double fi = -((double)p)*m_conus.dPhi;
			double degPhi = RadToDeg(-fi);
			outFile << endl << tt << " " << degPhi << " ";

			matrix m = M(p, t);
			outFile << m;
		}
	}
}

HandlerBackScatterPoint::HandlerBackScatterPoint(Particle *particle,
												 Light *incidentLight,
												 float wavelength)
	: HandlerPO(particle, incidentLight, wavelength)
{
}

void HandlerBackScatterPoint::HandleBeams(std::vector<Beam> &beams)
{
	Point3d vr(0, 0, 1);
	Point3d vf = -m_incidentLight->polarizationBasis;

#ifdef _DEBUG // DEB
	int c = 0;
#endif
	for (Beam &beam : beams)
	{
//		if (beam.direction.cz < BEAM_DIR_LIM)
//		{
//			continue;
//		}
#ifdef _DEBUG // DEB
		++c;
		vector<int> tr;
		Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
#endif
		int groupId = m_tracks->FindGroupByTrackId(beam.id);

		if (groupId < 0 && m_tracks->shouldComputeTracksOnly)
		{
			continue;
		}

		beam.polarizationBasis = beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);

		// absorbtion
		if (m_hasAbsorbtion && beam.nActs > 0)
		{
			ApplyAbsorbtion(beam);
		}

		Point3f beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
		beamBasis = beamBasis/Length(beamBasis);

		Point3f center = beam.Center();
		double projLenght = beam.opticalPath + DotProduct(center, beam.direction);

		matrixC jones(2, 2);
		matrixC fnJones = ComputeFnJones(beam.J, center, vr, projLenght);
		ApplyDiffraction(beam, beamBasis, vf, vr, fnJones, jones);

		// correction
		Matrix2x2c jonesCor = jones;
		jonesCor.m12 -= jonesCor.m21;
		jonesCor.m12 /= 2;
		jonesCor.m21 = -jonesCor.m12;

		if (groupId < 0 && !m_tracks->shouldComputeTracksOnly)
		{
			originContrib->AddToMueller(jones);
			correctedContrib->AddToMueller(jonesCor);
		}
		else
		{
			originContrib->AddToGroup(jones, groupId);
			correctedContrib->AddToGroup(jonesCor, groupId);
		}
	}

	originContrib->SumGroupTotal();
	correctedContrib->SumGroupTotal();
}


void HandlerBackScatterPoint::SetTracks(Tracks *tracks)
{
	Handler::SetTracks(tracks);
	originContrib = new PointContribution(tracks->size(), m_normIndex);
	correctedContrib = new PointContribution(tracks->size(), m_normIndex);
}

void HandlerBackScatterPoint::OutputContribution(ScatteringFiles &files,
												 double angle, double energy,
												 bool isOutputGroups,
												 string prefix)
{
	PointContribution *contrib = (prefix == "") ? originContrib : correctedContrib;
	contrib->SumTotal();

	energy *= m_normIndex;

	ofstream *all = files.GetMainFile(prefix + "all");
	*(all) << angle << ' ' << energy << ' ';
	*(all) << contrib->GetTotal() << endl;
//cout << endl << endl << contrib->GetRest()(0,0) << endl << endl ;
	if (isOutputGroups)
	{
		for (size_t gr = 0; gr < m_tracks->size(); ++gr)
		{
			ofstream &file = *(files.GetGroupFile(gr));
			file << angle << ' ' << energy << ' ';
			file << contrib->GetGroupMueller(gr) << endl;
		}
	}

	if (!m_tracks->shouldComputeTracksOnly)
	{
		ofstream &other = *(files.GetMainFile(prefix + "other"));
		other << angle << ' ' << energy << ' ';
		other << contrib->GetRest() << endl;

		ofstream &diff = *(files.GetMainFile(prefix + "difference"));
		diff << angle << ' ' << energy << ' ';
		diff << contrib->GetGroupTotal() << endl;
	}

	contrib->Reset();
}
