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
	  m_sphere(0.0, 0, 0),
	  m_nBadBeams(0)
{
	m_wavenumber = M_2PI/m_wavelength;
	m_wi2 = m_wavenumber*m_wavenumber;

	complex one(0, -1);
	m_complWave = (one * m_wavelength) / SQR(M_2PI);
	m_invComplWave = -one/m_wavelength;

	m_eps1 = 1e9*DBL_EPSILON;
	m_eps2 = 1e6*DBL_EPSILON;
	m_eps3 = 1e2;

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
		double d = m_scattering->MeasureOpticalPath(
					beam, beam.arr[i], tr);
		lengths.push_back(d);
	}

	Vector3f _n = beam.Normal();
	Point3d n = Point3d(_n.cx, _n.cy, _n.cz);

	Point3d hor;
	Point3d ver;
	ComputeCoordinateSystemAxes(n, hor, ver);

	Point3f cntr = beam.Center();
	Point3d center = ChangeCoordinateSystem(hor, ver, n,
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
	double path = m_scattering->MeasureOpticalPath(beam, beam.Center(), tr);

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

Point3d Handler::ChangeCoordinateSystem(const Point3d& hor, const Point3d& ver,
										const Point3d& normal,
										const Point3d& point) const
{
	// расчёт коор-т в СК наблюдателя
	const Point3d p_pr = point - normal*DotProductD(normal, point);
	return Point3d(DotProductD(p_pr, hor), DotProductD(p_pr, ver), 0);
}

void Handler::ComputeCoordinateSystemAxes(const Point3d& normal,
										  Point3d &hor, Point3d &ver) const
{
	if (fabs(normal.z) > 1-DBL_EPSILON)
	{
		hor = Point3d(0, -normal.z, 0);
		ver = Point3d(1, 0, 0);
	}
	else
	{
		const double tmp = sqrt(SQR(normal.x) + SQR(normal.y));
		hor = Point3d(normal.y/tmp, -normal.x/tmp, 0);
		ver = CrossProductD(normal, hor);
	}
}

void Handler::ComputeLengthIndices(const Beam &beam, BeamInfo &info)
{
	auto *lens = info.opticalLengths;

	Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.arr[0]) - info.projectedCenter;
	Point3d p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.arr[1]) - info.projectedCenter;
	Point3d p3 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.arr[2]) - info.projectedCenter;

	double den = p1.x*p2.y - p1.x*p3.y -
			p2.x*p1.y + p2.x*p3.y +
			p3.x*p1.y - p3.x*p2.y;

	if (fabs(den) < 1e-3)
	{
		m_isBadBeam = true;
		++m_nBadBeams;
	}

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

BeamInfo Handler::ComputeBeamInfo(const Beam &beam)
{
	BeamInfo info;
	info.normal = beam.Normal();
	info.normald = Point3d(info.normal.cx, info.normal.cy, info.normal.cz);

	bool isCcwOrder = DotProduct(info.normal, beam.direction) > 0;
	info.order.SetOrder(isCcwOrder, beam.nVertices);

	if (!isCcwOrder)
	{
		info.normal = -info.normal;
		info.normald = -info.normald;
	}

	ComputeCoordinateSystemAxes(info.normald, info.horAxis, info.verAxis);

	info.centerf = beam.Center();
	info.center = info.centerf;
	info.projectedCenter = ChangeCoordinateSystem(info.horAxis, info.verAxis,
												  info.normald, info.center);
	if (m_hasAbsorption && beam.lastFacetId != INT_MAX)
	{
		ComputeOpticalLengths(beam, info);
		ComputeLengthIndices(beam, info);
		m_cAbsExp = m_complWave * exp(m_cAbs*info.lenIndices.z);
	}

	info.area = beam.Area();

	info.projLenght = beam.opticalPath + DotProductD(info.center, beam.direction);

	info.beamBasis = CrossProduct(beam.polarizationBasis, beam.direction);
	info.beamBasis = info.beamBasis/Length(info.beamBasis); // basis of beam

	return info;
}

void Handler::ComputeOpticalLengths(const Beam &beam, BeamInfo &info)
{
	std::vector<int> tr;
	Tracks::RecoverTrack(beam, m_particle->nFacets, tr);
//	std::cout << tr.size() << std::endl;

	for (int i = 0; i < 3; ++i)
	{
		info.opticalLengths[i] = m_scattering->MeasureOpticalPath(
					beam, beam.arr[i], tr);
	}
	//	ExtropolateOpticalLenght(beam, tr);
}

Tracks *Handler::GetTracks() const
{
	return m_tracks;
}

Point3d Handler::ChangeCoordinateSystem(const Point3d& normal,
										const Point3d &point) const
{
	Point3d hor;  // условная горизонталь СК экрана в СК тела
	Point3d ver;  // третья ось (условная вертикаль СК экрана)
	ComputeCoordinateSystemAxes(normal, hor, ver);
	return ChangeCoordinateSystem(hor, ver, normal, point);
}

complex Handler::DiffractInclineAbs(const BeamInfo &info, const Beam &beam,
									const Point3d &direction) const
{
	const Point3f &dir = beam.direction;
	Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

	Point3d	pProj = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										   info.normald, k_k0);

	const complex A(pProj.x, info.lenIndices.x*m_riIm);
	const complex B(pProj.y, info.lenIndices.y*m_riIm);

	if (abs(A) < m_eps2 && abs(B) < m_eps2)
	{
		return m_invComplWave * info.area;
	}

	complex s(0, 0);

	const VertexOrder &order = info.order;

	Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.arr[order.begin]) - info.projectedCenter;
	Point3d p2;

	if (abs(B) > abs(A))
	{
		for (int i = order.startIndex; i != order.endIndex; i += order.inc)
		{
			p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.arr[i]) - info.projectedCenter;

			const double ci = (p1.y - p2.y)/(p1.x - p2.x);

			if (fabs(ci) > m_eps3)
			{
				p1 = p2;
				continue;
			}

			double bi = m_wavenumber*(p1.y - ci*p1.x);

			complex Ei = A + ci*B;
			complex tmp;

			if (abs(Ei) < m_eps1)
			{
				double mul = (p2.x*p2.x - p1.x*p1.x)/2.0;
				tmp = complex(-m_wi2*real(Ei)*mul,
							  m_wavenumber*(p2.x - p1.x) + m_wi2*imag(Ei)*mul);
			}
			else
			{
				double kReCi = m_wavenumber*real(Ei);
				double kImCi = -m_wavenumber*imag(Ei);
				tmp = (exp_im(kReCi*p2.x)*exp(kImCi*p2.x) -
					   exp_im(kReCi*p1.x)*exp(kImCi*p1.x))/Ei;
			}

			s += exp_im(bi*real(B)) * tmp * exp(-bi*imag(B));

#ifdef _DEBUG // DEB
		if (isnan(real(s)))
			int fff = 0;
#endif
			p1 = p2;
		}

		s /= B;
	}
	else
	{
		for (int i = order.startIndex; i != order.endIndex; i += order.inc)
		{
			p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.arr[i]) - info.projectedCenter;

			const double ci = (p1.x - p2.x)/(p1.y - p2.y);

			if (fabs(ci) > m_eps3)
			{
				p1 = p2;
				continue;
			}

			double bi = m_wavenumber*(p1.x - ci*p1.y);

			complex Ei = B + A*ci;
			complex tmp;

			if (abs(Ei) < m_eps1)
			{
				double mul = (p2.y*p2.y - p1.y*p1.y)/2.0;
				tmp = complex(-m_wi2*real(Ei)*mul,
							  m_wavenumber*(p2.y - p1.y) + m_wi2*imag(Ei)*mul);
			}
			else
			{
				double kReEi = m_wavenumber*real(Ei);
				double kImEi = -m_wavenumber*imag(Ei);
				tmp = (exp_im(kReEi*p2.y)*exp(kImEi*p2.y) -
					   exp_im(kReEi*p1.y)*exp(kImEi*p1.y))/Ei;
			}

			s += exp_im(bi*real(A)) * exp(-bi*imag(A)) * tmp;

#ifdef _DEBUG // DEB
		if (isnan(real(s)))
			int fff = 0;
#endif
			p1 = p2;
		}

		s /= -A;
	}

#ifndef _DEBUG // DEB
//	double dddd = exp(m_absMag*info.lenIndices.z);
#endif
	return s * m_cAbsExp;
}

double Handler::BeamCrossSection(const Beam &beam) const
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

complex Handler::DiffractIncline(const BeamInfo &info, const Beam &beam,
								 const Point3d &direction) const
{
	const Point3f &dir = beam.direction;
	Point3d k_k0 = -direction + Point3d(dir.cx, dir.cy, dir.cz);

	Point3d	pt_proj = ChangeCoordinateSystem(info.horAxis, info.verAxis,
											 info.normald, k_k0);
	const double A = pt_proj.x;
	const double B = pt_proj.y;

	double absA = fabs(A);
	double absB = fabs(B);

	if (absA < m_eps2 && absB < m_eps2)
	{
		return m_invComplWave * info.area;
	}

	complex s(0, 0);
	const VertexOrder &order = info.order;

	Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.arr[order.begin]) - info.projectedCenter;
	Point3d p2;

	if (absB > absA)
	{
		for (int i = order.startIndex; i != order.endIndex; i += order.inc)
		{
			p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.arr[i]) - info.projectedCenter;

			if (fabs(p1.x - p2.x) < m_eps1)
			{
				p1 = p2;
				continue;
			}

			const double ai = (p1.y - p2.y)/(p1.x - p2.x);
			const double Ci = A+ai*B;

			complex tmp;

			if (fabs(Ci) < m_eps1)
			{
				tmp = complex(-m_wi2*Ci*(p2.x*p2.x - p1.x*p1.x)/2.0,
							  m_wavenumber*(p2.x - p1.x));
			}
			else
			{
				double kCi = m_wavenumber*Ci;
				tmp = (exp_im(kCi*p2.x) - exp_im(kCi*p1.x))/Ci;
			}

			const double bi = p1.y - ai*p1.x;
			s += exp_im(m_wavenumber*B*bi) * tmp;

			p1 = p2;
		}

		s /= B;
	}
	else
	{
		for (int i = order.startIndex; i != order.endIndex; i += order.inc)
		{
			p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.arr[i]) - info.projectedCenter;

			if (fabs(p1.y - p2.y) < m_eps1)
			{
				p1 = p2;
				continue;
			}

			const double ci = (p1.x - p2.x)/(p1.y - p2.y);
			const double Ei = A*ci+B;

			complex tmp;

			if (fabs(Ei) < m_eps1)
			{
				tmp = complex(-m_wi2*Ei*(p2.y*p2.y - p1.y*p1.y)/2.0,
							  m_wavenumber*(p2.y - p1.y));
			}
			else
			{
				double kEi = m_wavenumber*Ei;
				tmp = (exp_im(kEi*p2.y) - exp_im(kEi*p1.y))/Ei;
			}

			const double di = p1.x - ci*p1.y;
			s += exp_im(m_wavenumber*A*di) * tmp;

			p1 = p2;
		}

		s /= -A;
	}

	return m_complWave * s;
}
