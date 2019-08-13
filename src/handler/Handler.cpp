#include "Handler.h"

#include "Mueller.hpp"
#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

Handler::Handler(Particle *particle, Light *incidentLight, double wavelength)
	: m_incidentLight(incidentLight),
	  m_particle(particle),
	  m_wavelength(wavelength),
	  m_hasAbsorption(false),
	  m_normIndex(1),
	  m_sphere(0.0, 0, 0)
{
//	m_wavelength = 0.532;
	m_wavenumber = M_2PI/m_wavelength;
	m_wn2 = m_wavenumber*m_wavenumber;

	complex one(0, -1);
	m_complWave = (one * m_wavelength) / SQR(M_2PI);
	m_invComplWave = -one/m_wavelength;

	m_eps1 = 1e9*DBL_EPSILON;
	m_eps2 = 1e6*DBL_EPSILON;
	m_eps3 = 1e-2;
//	m_eps3 = m_eps1;


	m_logFile.open("log1.txt", std::ios::out);
	m_logFile << std::setprecision(std::numeric_limits<long double>::digits10 + 1);
}

void Handler::HandleBeams(std::vector<Beam> &/*beams*/)
{
}

void Handler::SetTracks(Tracks *tracks)
{
//	if (!m_tracks)
//	{
//		std::cerr << "Tracks are not set" << std::endl;
//		throw std::exception();
//	}

	m_tracks = tracks;
}

void Handler::WriteMatricesToFile(string &/*destName*/)
{
}

void Handler::SetNormIndex(double value)
{
	m_normIndex = value;
}

void Handler::SetSinZenith(double value)
{
	m_sinZenith = value;
}

void Handler::SetAbsorptionAccounting(bool value)
{
	m_hasAbsorption = value;
	m_ri = m_particle->GetRefractiveIndex();
	m_riIm = imag(m_ri);
	m_cAbs = -m_riIm*m_wavenumber;
}

void Handler::SetScattering(Scattering *scattering)
{
	m_scattering = scattering;
}

void Handler::ExtropolateOpticalLenght(Beam &beam, const std::vector<int> &tr)
{
	std::vector<double> lengths;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		double d = m_scattering->ComputeInternalOpticalPath(
					beam, beam.arr[i], tr);
		lengths.push_back(d);
	}

	Vector3f _n = beam.Normal();
	Point3d n = Point3d(_n.cx, _n.cy, _n.cz);

	Axes csAxes;
	ComputeCoordinateSystemAxes(n, csAxes);

	Point3f cntr = beam.Center();
	Point3d center = ChangeCoordinateSystem(csAxes, n,
											Point3d(cntr.cx, cntr.cy, cntr.cz));
	double ls[3];
	ls[0] = lengths[0];
	ls[1] = lengths[1];
	ls[2] = lengths[2];

//	Point3d lens = ComputeLengthIndices(beam, info);

//	for (int i = 3; i < beam.nVertices; ++i)
//	{
//		Point3d newP = ChangeCoordinateSystem(hor, ver, n, beam.arr[i]) - center;
//		double newL = lens.z + lens.x*newP.x + lens.y*newP.y;
//		double err = fabs((lengths[i] - newL)/lengths[i])*100;
//		if (err > 5)
//			m_logFile << Polygon(beam) << "Area: " << beam.Area() << ' '
//					  << i << ", " << "Error: " << err << std::endl;
//	}
}

void Handler::ApplyAbsorption(Beam &beam)
{
	vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);

//	double opAbs = CalcOpticalPathAbsorption(beam);
	double path = m_scattering->ComputeInternalOpticalPath(beam, beam.Center(), tr);

#ifdef _DEBUG // DEB
//	double ddd = fabs(path - beam.opticalPath);
//	if (fabs(path - beam.opticalPath) >= 10e-4)
//		int ggg = 0;
#endif

	if (path > DBL_EPSILON)
	{
		double abs = exp(m_cAbs*path);
		beam.J *= abs;
	}
}

Point3d Handler::ChangeCoordinateSystem(const Axes &axes, const Point3d& normal,
										const Point3d& point) const
{
	// расчёт коор-т в СК наблюдателя
	const Point3d p_pr = point - normal*DotProductD(normal, point);

	return Point3d(DotProductD(p_pr, axes.horisontal),
				   DotProductD(p_pr, axes.vertical), 0);
}

void Handler::ComputeCoordinateSystemAxes(const Point3d& normal,
										  Axes &axes) const
{
	if (fabs(normal.z) > 1-DBL_EPSILON)
	{
		axes.horisontal = Point3d(0, -normal.z, 0);
		axes.vertical = Point3d(1, 0, 0);
	}
	else
	{
		const double tmp = sqrt(SQR(normal.x) + SQR(normal.y));
		axes.horisontal = Point3d(normal.y/tmp, -normal.x/tmp, 0);
		axes.vertical = CrossProductD(normal, axes.horisontal);
	}
}

void Handler::ComputeLengthIndices(const Beam &beam, BeamInfo &info) const
{
	auto *lens = info.opticalLengths;

	Point3d p1 = ChangeCoordinateSystem(info.csAxes, info.normald, beam.arr[0])
			- info.projectedCenter;
	Point3d p2 = ChangeCoordinateSystem(info.csAxes, info.normald, beam.arr[1])
			- info.projectedCenter;
	Point3d p3 = ChangeCoordinateSystem(info.csAxes, info.normald, beam.arr[2])
			- info.projectedCenter;

	double den = p1.x*p2.y - p1.x*p3.y -
			p2.x*p1.y + p2.x*p3.y +
			p3.x*p1.y - p3.x*p2.y;

	info.lenIndices.z = (lens[0]*p2.x*p3.y - lens[0]*p3.x*p2.y -
			lens[1]*p1.x*p3.y + lens[1]*p3.x*p1.y +
			lens[2]*p1.x*p2.y - lens[2]*p2.x*p1.y) / den;

	info.lenIndices.x = (lens[0]*p2.y - lens[0]*p3.y -
			lens[1]*p1.y + lens[1]*p3.y +
			lens[2]*p1.y - lens[2]*p2.y) / den;

	info.lenIndices.y = -(lens[0]*p2.x - lens[0]*p3.x -
			lens[1]*p1.x + lens[1]*p3.x +
			lens[2]*p1.x - lens[2]*p2.x) / den;
}

Tracks *Handler::GetTracks() const
{
	return m_tracks;
}

Point3d Handler::ChangeCoordinateSystem(const Point3d& normal,
										const Point3d &point) const
{
	Axes csAxes;  /* условная горизонталь СК экрана в СК тела и
					третья ось (условная вертикаль СК экрана)*/
	ComputeCoordinateSystemAxes(normal, csAxes);
	return ChangeCoordinateSystem(csAxes, normal, point);
}

complex Handler::ComputeAvgBeamEnergy(const Polygon &pol, const BeamInfo &info,
									  Point3d &p1, Point3d &p2,
									  double &p11, double &p12,
									  double &p21, double &p22,
									  const complex &c1, const complex &c2) const
{
	complex s(0, 0);
	complex tmp;
	double ai;
	complex ci;

	for (int i = info.order.startIndex; i != info.order.endIndex; i += info.order.increment)
	{
		p2 = ChangeCoordinateSystem(info.csAxes, info.normald, pol.arr[i])
				- info.projectedCenter;

		if (fabs(p11 - p21) > m_eps3)
		{
			ai = (p12 - p22)/(p11 - p21);
			ci = c1 + ai*c2;

			if (abs(ci) < m_eps1)
			{
				double mul = (p21*p21 - p11*p11)/2.0;
				tmp = mul*complex(-m_wn2*real(ci),
								  m_wavenumber*(p21 - p11) + m_wn2*imag(ci));
			}
			else
			{
				double cRe = m_wavenumber*real(ci);
				double cIm = -m_wavenumber*imag(ci);
				tmp = (exp_im(cRe*p21)*exp(cIm*p21) -
					   exp_im(cRe*p11)*exp(cIm*p11))/ci;
			}

			double kDi = m_wavenumber*(p12 - ai*p11);
			s += exp_im(kDi*real(c2)) * exp(-kDi*imag(c2)) * tmp;
		}

		p1 = p2;
	}

	s /= c2;
	return s;
}

complex Handler::DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
									const Point3d &direction) const
{
	const Point3f &beamDir = beam.direction;
	Point3d k_k0 = -direction + Point3d(beamDir.cx, beamDir.cy, beamDir.cz);

	Point3d	pProj = ChangeCoordinateSystem(info.csAxes, info.normald, k_k0);

	const complex A(pProj.x, info.lenIndices.x*m_riIm);
	const complex B(pProj.y, info.lenIndices.y*m_riIm);

	if (abs(A) > m_eps2 || abs(B) > m_eps2)
	{
		Point3d p1 = ChangeCoordinateSystem(info.csAxes, info.normald,
											beam.arr[info.order.begin])
				- info.projectedCenter;
		Point3d p2;

		complex s(0, 0);

		if (abs(B) > abs(A))
		{
			s = ComputeAvgBeamEnergy(beam, info, p1, p2, p1.x, p2.x, p1.y, p2.y,
									 A, B);
		}
		else
		{
			s = -ComputeAvgBeamEnergy(beam, info, p1, p2, p1.y, p2.y, p1.x, p2.x,
									  B, A);
		}

#ifdef _DEBUG // DEB
		double dddd = exp(m_cAbs*info.lenIndices.z);
#endif
		return m_complWave * s * exp(m_cAbs*info.lenIndices.z);
	}
	else
	{
		return m_invComplWave * info.area;
	}
}

double Handler::BeamCrossSection(const Beam &beam) const
{
	Point3f normal = m_particle->facets[beam.lastFacetId].ex_normal; // normal of last facet of beam
	double cosA = DotProduct(normal, beam.direction);
	double e = fabs(cosA);

	if (e > m_eps1)
	{
		double area = beam.Area();
		double len = Length(normal);
		return (e*area)/len;
	}
	else
	{
		return 0;
	}
}

complex Handler::DiffractIncline(const BeamInfo &info, const Beam &beam,
								 const Point3d &direction) const
{
	const Point3f &dir = beam.direction;
	Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

	Point3d	pt_proj = ChangeCoordinateSystem(info.csAxes, info.normald, k_k0);
	const double A = pt_proj.x;
	const double B = pt_proj.y;

	double absA = fabs(A);
	double absB = fabs(B);

	if (absA < m_eps2 && absB < m_eps2)
	{
		return m_invComplWave * info.area;
	}

	complex s(0, 0);

	Point3d p1 = ChangeCoordinateSystem(info.csAxes, info.normald,
										beam.arr[info.order.begin])
			- info.projectedCenter;
	Point3d p2;

	if (absB > absA)
	{
		for (int i = info.order.startIndex; i != info.order.endIndex; i += info.order.increment)
		{
			p2 = ChangeCoordinateSystem(info.csAxes, info.normald, beam.arr[i])
					- info.projectedCenter;

			if (fabs(p1.x - p2.x) > m_eps1)
			{
				const double ai = (p1.y - p2.y)/(p1.x - p2.x);
				const double Ci = A+ai*B;

				complex tmp;

				if (fabs(Ci) < m_eps1)
				{
					tmp = complex(-m_wn2*Ci*(p2.x*p2.x - p1.x*p1.x)/2.0,
								  m_wavenumber*(p2.x - p1.x));
				}
				else
				{
					double kCi = m_wavenumber*Ci;
					tmp = (exp_im(kCi*p2.x) - exp_im(kCi*p1.x))/Ci;
				}

				const double bi = p1.y - ai*p1.x;
				s += exp_im(m_wavenumber*B*bi) * tmp;
			}

			p1 = p2;
		}

		s /= B;
	}
	else
	{
		for (int i = info.order.startIndex; i != info.order.endIndex; i += info.order.increment)
		{
			p2 = ChangeCoordinateSystem(info.csAxes, info.normald, beam.arr[i])
					- info.projectedCenter;

			if (fabs(p1.y - p2.y) > m_eps1)
			{
				const double ci = (p1.x - p2.x)/(p1.y - p2.y);
				const double Ei = A*ci+B;

				complex tmp;

				if (fabs(Ei) < m_eps1)
				{
					tmp = complex(-m_wn2*Ei*(p2.y*p2.y - p1.y*p1.y)/2.0,
								  m_wavenumber*(p2.y - p1.y));
				}
				else
				{
					double kEi = m_wavenumber*Ei;
					tmp = (exp_im(kEi*p2.y) - exp_im(kEi*p1.y))/Ei;
				}

				const double di = p1.x - ci*p1.y;
				s += exp_im(m_wavenumber*A*di) * tmp;
			}

			p1 = p2;
		}

		s /= -A;
	}

	return m_complWave * s;
}
