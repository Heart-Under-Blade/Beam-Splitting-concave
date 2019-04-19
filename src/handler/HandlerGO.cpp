#include "HandlerGO.h"

#include "Particle.h"
#include "Mueller.hpp"

#include <limits>
#include <iomanip>

#define BIN_SIZE			M_PI/SPHERE_RING_NUM

using namespace std;

HandlerGO::HandlerGO(Particle *particle, Light *incidentLight, float wavelength)
	: Handler(particle, incidentLight, wavelength)
{
	m_logFile.open("log1.txt", ios::out);
	m_logFile << setprecision(numeric_limits<long double>::digits10 + 1);
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

void HandlerGO::WriteLog(const string &str)
{
	m_logFile << str;
}

void HandlerGO::MultiplyMueller(const Beam &beam, matrix &m)
{
	double cross = BeamCrossSection(beam);
	double area = cross*m_sinAngle;
	m *= area;
}

void HandlerGO::WriteMatricesToFile(string &destName)
{
	AverageOverAlpha(true, m_normIndex, m_totalContrib, destName);
	WriteToFile(m_totalContrib, m_normIndex, destName + "_all");
}

matrix HandlerGO::ComputeMueller(int zenAng, Beam &beam)
{
	matrix m = Mueller(beam.Jones);

	if (zenAng < 180 && zenAng > 0)
	{
		const float &y = beam.direction.coordinates[1];

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
	const float &x = dir.coordinates[0];
	const float &y = dir.coordinates[1];

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

	allFile << "tetta"\
				   "M11 M12/M11 M13/M11 M14/M11 "\
			   "M21/M11 M22/M11 M23/M11 M24/M11 "\
			   "M31/M11 M32/M11 M33/M11 M34/M11 "\
			   "M21/M41 M42/M11 M43/M11 M44/M11 ";

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

Point3f HandlerGO::CalcK(vector<int> &tr)
{	// OPT: сделать из переменных ссылки
	Point3f k, tmp;
	Point3f n1 = m_particle->GetActualFacet(tr[0])->in_normal;
	Point3f nq = m_particle->GetActualFacet(tr[tr.size()-1])->in_normal;
	Point3f::CrossProduct(nq, n1, tmp);
	Point3f::CrossProduct(tmp, nq, k);

	for (int i = tr.size()-2; i > 0; --i)
	{
		Point3f ni = m_particle->GetActualFacet(tr[i])->in_normal;
		k = k - ni*2*fabs(Point3f::DotProduct(ni, k));
	}

	Point3f::Normalize(k);
	return k;
}

double HandlerGO::ComputeOpticalPathAbsorption(const Beam &beam)
{	// OPT: вынести переменные из цикла
	double opticalPath = 0;

	vector<int> tr;
	m_tracks->RecoverTrack(beam, tr);

	Point3f k = CalcK(tr);

	Point3f n1 = m_particle->GetActualFacet(tr[0])->in_normal;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		double delta = Point3f::Length(beam.Center() - beam.vertices[i])/Point3f::Length(k);
		opticalPath += (delta*Point3f::DotProduct(k, n1))/Point3f::DotProduct(beam.direction, n1);
	}

	opticalPath /= beam.nVertices;

	return opticalPath;
}
