#include "Handler.h"
#include "Mueller.hpp"
#include <iostream>
#include <limits>
#include <iomanip>

#define SPHERE_RING_NUM		180		// number of rings no the scattering sphere
#define BIN_SIZE			M_PI/SPHERE_RING_NUM

using namespace std;

HandlerGO::HandlerGO(Particle *particle, Light *incidentLight, float wavelength)
	: m_particle(particle),
	  m_incidentLight(incidentLight),
	  m_wavelength(wavelength),
	  m_isAccountAbsorbtion(false)
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

	f[0] =  contrib.forw[0][0];
	f[1] = (contrib.forw[1][1] + contrib.forw[2][2])/2.0;
	f[2] =  contrib.forw[3][3];

	// Extracting the forward and backward peak in a separate file if needed
	if (EDF)
	{
		ExtractPeaks(b, f, norm);
	}
	else
	{
		contrib.scatMatrix(0,SPHERE_RING_NUM,0,0) += f[0];
		contrib.scatMatrix(0,0,0,0) += b[0];
		contrib.scatMatrix(0,SPHERE_RING_NUM,1,1) += f[1];
		contrib.scatMatrix(0,0,1,1) += b[1];
		contrib.scatMatrix(0,SPHERE_RING_NUM,2,2) += f[1];
		contrib.scatMatrix(0,0,2,2) -= b[1];
		contrib.scatMatrix(0,SPHERE_RING_NUM,3,3) += f[2];
		contrib.scatMatrix(0,0,3,3) += b[2];
	}
}

void HandlerGO::HandleBeams(std::vector<Beam> &beams, double angle)
{
}

double HandlerGO::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f normal = m_particle->facets[beam.lastFacetID].ex_normal; // normal of last facet of beam
	double cosFB = DotProduct(normal, beam.light.direction);
	double e = fabs(cosFB);

	if (e < eps)
	{
		return 0;
	}

	double area = beam.Area();
	double n = Length(normal);
	return (e*area)/n;
}

void HandlerGO::AddToResultMullerGO(const Point3f &dir, matrix &bf, double area,
									ContributionGO &contr)
{
	const float &z = dir.cz;

	// Collect the beam in array
	if (z >= 1-DBL_EPSILON)
	{
		contr.back += area*bf;
	}
	else if (z <= DBL_EPSILON-1)
	{
		contr.forw += area*bf;
	}
	else
	{
		const float &y = dir.cy;

		if (y*y > DBL_EPSILON)
		{	// rotate the Mueller matrix of the beam to appropriate coordinate system
			RotateMuller(dir, bf);
		}

#ifdef _CALC_AREA_CONTIBUTION_ONLY
		bf = matrix(4,4);
		bf.Identity();
#endif
		const unsigned int zenAng = round(acos(z)/(M_PI/SPHERE_RING_NUM));
		contr.scatMatrix.insert(0, zenAng, area*bf);
	}
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

		matrix bf = contrib.scatMatrix(0, j);

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

double HandlerGO::CalcOpticalPathAbsorption(const Beam &beam)
{	// OPT: вынести переменные из цикла
	double opticalPath = 0;

	vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->facetNum, tr);

	Point3f k = CalcK(tr);

//	if (Length(k) > FLT_EPSILON)
//	{
		Point3f n1 = m_particle->facets[tr[0]].in_normal;

		for (int i = 0; i < beam.size; ++i)
		{
//			double lk = Length(k);
			double delta = Length(beam.Center() - beam.arr[i])/*/lk*/;
			double a1 = (delta*DotProduct(k, n1));
			double a2 = DotProduct(beam.light.direction, n1);
			opticalPath += a1/a2;
		}

		opticalPath /= beam.size;
//	}

	return opticalPath;
}

void HandlerGO::WriteMatricesToFile(double norm, string &filename)
{
}

void HandlerGO::SetTracks(Tracks *tracks)
{
	if (!m_tracks)
	{
		std::cerr << "Tracks are not set" << std::endl;
		throw std::exception();
	}

	m_tracks = tracks;
	m_groupContrib.resize(m_tracks->size());
}

double HandlerGO::CalcTotalScatteringEnergy(double norm)
{
	double D_tot = m_totalContrib.back[0][0] + m_totalContrib.forw[0][0];

	for (int i = 0; i <= SPHERE_RING_NUM; ++i)
	{
		D_tot += m_totalContrib.scatMatrix(0, i, 0, 0);
	}

	return D_tot*norm;
}

void HandlerGO::SetAbsorbtionAccounting(bool value)
{
	m_isAccountAbsorbtion = value;
}

void HandlerGO::WriteLog(const string &str)
{
	m_logFile << str;
}

void HandlerGO::SetTracing(Tracing *tracing)
{
	m_tracing = tracing;
}


HandlerTotalGO::HandlerTotalGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerTotalGO::HandleBeams(std::vector<Beam> &beams, double angle)
{
	double sinBeta = sin(angle);
#ifdef _DEBUG // DEB
	int c = 0;
#endif
	double kAbs = -M_2PI*imag(m_particle->GetRefractiveIndex())/m_wavelength;

	for (Beam &beam : beams)
	{
		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);

		double cross = BeamCrossSection(beam);
		double area = cross*sinBeta;

		// absorbtion
		if (/*m_isAccountAbsorbtion &&*/ beam.level > 0)
		{
			vector<int> tr;
			Tracks::RecoverTrack(beam, m_particle->facetNum, tr);

//			double opAbs = CalcOpticalPathAbsorption(beam);
			double path = m_tracing->ComputeInternalOpticalPath(beam, tr);
#ifdef _DEBUG // DEB
if (/*beam.level > 1 && */isnan(path))
{
	double s = 0;
}

m_logFile << fabs(path - beam.opticalPath) << endl;
//m_logFile << path << ';' << beam.opticalPath << endl;
//cout << endl << endl << "path " << path << endl << "opticalPath " << beam.opticalPath;
#endif
			if (fabs(path) > DBL_EPSILON)
			{
				double abs = exp(kAbs*path);
				beam.J *= abs;
			}
		}

		matrix bf = Mueller(beam.J);
#ifdef _DEBUG // DEB
if (beam.level > 2)
{
	double s = bf[0][0];
//	double f = Mueller(Jd)[0][0];
//	m_logFile << s << endl;
}
#endif
//		if (isinf(bf[0][0]))
//			int f = 0;
		AddToResultMullerGO(beam.light.direction, bf, area, m_totalContrib);
#ifdef _DEBUG // DEB
//		double oo =  m_totalContrib.scatMatrix(0, 6)[0][0];
//		m_logFile << c++ << " " << oo << endl;
#endif
	}

// deb
//m_logFile << m_totalContrib.scatMatrix(0, 6)[0][0] << endl;

	beams.clear();
}

HandlerGroupGO::HandlerGroupGO(Particle *particle, Light *incidentLight, float wavelength)
	: HandlerGO(particle, incidentLight, wavelength)
{
}

void HandlerGroupGO::HandleBeams(std::vector<Beam> &beams, double angle)
{
	double sinBeta = sin(angle);

	for (Beam &beam : beams)
	{
		int grID = m_tracks->FindGroup(beam.id);

		if (grID < 0)
		{
			continue;
		}

		beam.RotateSpherical(-m_incidentLight->direction,
							 m_incidentLight->polarizationBasis);

		double cross = BeamCrossSection(beam);
		double area = cross*sinBeta;
		matrix bf = Mueller(beam.J);

		// group contribution
		AddToResultMullerGO(beam.light.direction, bf, area, m_groupContrib[grID]);

		// total contribution
		AddToResultMullerGO(beam.light.direction, bf, area, m_totalContrib);
	}

	beams.clear();
}

void HandlerGroupGO::WriteMatricesToFile(double norm, string &dirName)
{
	string dir = CreateFolder(dirName);

	for (size_t i = 0; i < m_groupContrib.size(); ++i)
	{
		if ((*m_tracks)[i].size != 0)
		{
			string subname = (*m_tracks)[i].CreateGroupName();
			AverageOverAlpha(0, norm, m_groupContrib[i]);
			WriteToFile(m_groupContrib[i], norm, dir + subname);
		}
	}

	WriteMatricesToFile(norm, dirName);
}


void HandlerTotalGO::WriteMatricesToFile(double norm, string &filename)
{
	filename += "_all";
	AverageOverAlpha(0, norm, m_totalContrib);
	WriteToFile(m_totalContrib, norm, filename);
}
