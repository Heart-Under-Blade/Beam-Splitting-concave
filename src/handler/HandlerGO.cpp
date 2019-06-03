#include "HandlerGO.h"

#include <limits>
#include <iostream>
#include <iomanip>

#include "Mueller.hpp"

HandlerGO::HandlerGO(Particle *particle, Light *incidentLight, double wavelength)
	: Handler(particle, incidentLight, wavelength)
{
}

void HandlerGO::SetTracks(Tracks *tracks)
{
	Handler::SetTracks(tracks);
	m_tracksContrib.resize(m_tracks->size());
}

void HandlerGO::ExtractPeaks(double *b, double *f, double norm,
							 const std::string &destDir)
{
	std::ofstream bck(destDir + "_back.dat", std::ios::out);
	std::ofstream frw(destDir + "_forward.dat", std::ios::out);
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

void HandlerGO::AverageOverAlpha(int EDF, double norm, ContributionGO &contrib,
								 const std::string &destDir)
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
		ExtractPeaks(b, f, norm, destDir);
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

void HandlerGO::MultiplyMueller(const Beam &beam, matrix &m)
{
	double cross = BeamCrossSection(beam);
	double area = cross*m_sinZenith;
	m *= area;
}

matrix HandlerGO::ComputeMueller(float zenAng, Beam &beam)
{
	matrix m = Mueller(beam.J);
#ifdef _DEBUG // DEB
	double &ddd = m[0][0];
#endif
	if (zenAng < 180-FLT_EPSILON && zenAng > FLT_EPSILON)
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
	std::string name = CreateUniqueFileName(filename);
	std::ofstream allFile(name, std::ios::out);

	allFile << "tetta dS M11 M12/M11 M13/M11 M14/M11 "\
				"M21/M11 M22/M11 M23/M11 M24/M11 "\
				"M31/M11 M32/M11 M33/M11 M34/M11 "\
				"M41/M11 M42/M11 M43/M11 M44/M11";

	for (int j = SPHERE_RING_NUM; j >= 0; j--)
	{
		double tmp0 = 180.0/SPHERE_RING_NUM*(SPHERE_RING_NUM-j);
		double tmp1 = (j == 0) ? -(0.25*180.0)/SPHERE_RING_NUM : 0;
		double tmp2 = (j == (int)SPHERE_RING_NUM) ? (0.25*180.0)/SPHERE_RING_NUM : 0;

		double sn = (j == 0 || j == (int)SPHERE_RING_NUM)
				? 1-cos(BIN_SIZE/2.0)
				: (cos((j-0.5)*BIN_SIZE)-cos((j+0.5)*BIN_SIZE));

		// Special case in first and last step
		allFile << '\n' << tmp0 + tmp1 + tmp2 << ' ' << (M_2PI*sn);

		matrix bf = contrib.muellers(0, j);

#ifdef _DEBUG // DEB
		double dd = bf[0][0];
		std::cout << sn;
#endif

		if (bf[0][0] <= DBL_EPSILON)
		{
			allFile << " 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0";
		}
		else
		{
			allFile << ' ' << bf[0][0]*norm/(M_2PI*sn)
					<< ' ' << bf[0][1]/bf[0][0]
					<< ' ' << bf[0][2]/bf[0][0] // usually = 0
					<< ' ' << bf[0][3]/bf[0][0] // usually = 0
					<< ' ' << bf[1][0]/bf[0][0]
					<< ' ' << bf[1][1]/bf[0][0]
					<< ' ' << bf[1][2]/bf[0][0] // usually = 0
					<< ' ' << bf[1][3]/bf[0][0] // usually = 0
					<< ' ' << bf[2][0]/bf[0][0] // usually = 0
					<< ' ' << bf[2][1]/bf[0][0] // usually = 0
					<< ' ' << bf[2][2]/bf[0][0]
					<< ' ' << bf[2][3]/bf[0][0]
					<< ' ' << bf[3][0]/bf[0][0] // usually = 0
					<< ' ' << bf[3][1]/bf[0][0] // usually = 0
					<< ' ' << bf[3][2]/bf[0][0]
					<< ' ' << bf[3][3]/bf[0][0];
		}
	}

	allFile.close();
}

Point3f HandlerGO::CalcK(std::vector<int> &tr)
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

	std::vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

	Point3f k = CalcK(tr);
	Point3f n1 = m_particle->facets[tr[0]].in_normal;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		double delta = Length(beam.Center() - beam.arr[i])/Length(k);
		opticalPath += (delta*DotProduct(k, n1))/DotProduct(beam.direction, n1);
	}

	opticalPath /= beam.nVertices;
	return opticalPath;
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

void HandlerGO::WriteLog(const std::string &str)
{
	m_logFile << str;
}
