#include "Handler.h"

#include <iostream>
#include <limits>
#include <iomanip>
#include <algorithm>

#include "Mueller.hpp"

using namespace std;

Handler::Handler(Particle *particle, Light *incidentLight, double wavelength)
	: m_incidentLight(incidentLight),
	  m_particle(particle),
	  m_wavelength(wavelength),
	  m_hasAbsorption(false),
	  m_normIndex(1)
{
	m_wavenumber = M_2PI/m_wavelength;
	m_wi2 = m_wavenumber*m_wavenumber;

	complex one(0, -1);
	m_complWave = (one * m_wavelength) / SQR(M_2PI);
	m_invComplWave = -one/m_wavelength;

	m_eps1 = 1e9*DBL_EPSILON;
	m_eps2 = 1e6*DBL_EPSILON;

	m_riIm = imag(m_ri);
	m_absMag = -m_wavenumber*m_riIm;

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

void Handler::ExtropolateOpticalLenght(Beam &beam, const std::vector<int> &tr)
{
	std::vector<double> lengths;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		OpticalPath op = m_scattering->ComputeOpticalPath
				(beam, beam.vertices[i], tr);
		lengths.push_back(op.GetTotal());
	}

	Vector3f _n = beam.Normal();
	Point3d n = Point3d(_n.coordinates[0], _n.coordinates[1], _n.coordinates[2]);

	Point3d hor;
	Point3d ver;
	ComputeCoordinateSystemAxes(n, hor, ver);

	Point3f cntr = beam.Center();
	Point3d center = ChangeCoordinateSystem
			(hor, ver, n, Point3d(cntr.coordinates[0], cntr.coordinates[1],
			cntr.coordinates[2]));
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

double Handler::BeamCrossSection(const Beam &beam) const
{
	const double eps = 1e7*DBL_EPSILON;

	Point3f &normal = beam.facet->ex_normal; // normal of last facet of beam
	double cosA = Point3f::DotProduct(normal, beam.direction);
	double e = fabs(cosA);

	if (e < eps)
	{
		return 0;
	}

	double area = beam.Area();
	double len = Point3f::Length(normal);
	return (e*area)/len;
}

void Handler::WriteMatricesToFile(string &/*destName*/)
{
}

void Handler::SetNormIndex(double normIndex)
{
	m_normIndex = normIndex;
}

void Handler::SetSinZenith(double value)
{
	m_sinZenith = value;
}

void Handler::EnableAbsorption(bool isNew)
{
	m_hasAbsorption = true;
	m_isNewAbs = isNew;
	m_cAbs = -m_wavenumber*imag(m_scattering->m_refractiveIndex);

	m_absLogFile.open("abslog1.txt", ios::out);
	m_absLogFile << setprecision(10);
	m_absLogFile << "No" << ' '
				 << "CtrPath" << ' '
				 << "AvgPath" << ' '
				 << "Max-Min" << ' '
				 << "Nact" << ' '
				 << "Npt" << ' '
				 << "Tr" << ' '
				 << endl;
}

void Handler::SetScattering(Scattering *scattering)
{
	m_scattering = scattering;
}

void Handler::OutputPaths(const Beam &beam, const OpticalPath &path)
{
	vector<int> track;
	m_tracks->RecoverTrack(beam, track);

	vector<double> ps;
	double sum = 0;

	for (int i = 0; i < beam.nVertices; ++i)
	{
		OpticalPath p0 = m_scattering->ComputeOpticalPath(beam, beam.vertices[i],
														  track);
		ps.push_back(p0.internal);
		sum += p0.internal;
	}

	double maxPath = *std::max_element(ps.begin(), ps.end());
	double minPath = *std::min_element(ps.begin(), ps.end());
	double delta = fabs(path.internal - sum/ps.size());

	if (delta >= 10e-4)
	{
		m_absLogFile << ++count << ' '
					 << path.internal << ' '
					 << sum/ps.size() << ' '
					 << maxPath - minPath << ' '
					 << beam.actNo << ' '
					 << beam.nVertices << ' '
					 << Tracks::TrackToStr(track)
					 << endl;
	}
}

OpticalPath Handler::ComputeOpticalPath(const Beam &beam)
{
	vector<int> track;
	m_tracks->RecoverTrack(beam, track);
//	double opAbs = CalcOpticalPathAbsorption(beam);
	return m_scattering->ComputeOpticalPath(beam, beam.Center(), track);
}

void Handler::ApplyAbsorption(Beam &beam)
{
	auto path = ComputeOpticalPath(beam);

	if (path.internal > DBL_EPSILON)
	{
#ifdef _DEBUG // DEB
//		OutputPaths(beam, path);
//		if (fabs(path.GetTotal() - beam.opticalPath) >= 10e-4)
//			int ggg = 0;
#endif
		double abs = exp(m_cAbs*path.internal);
		beam.Jones *= abs;
	}
}

Point3d Handler::ChangeCoordinateSystem(const Point3d& hor, const Point3d& ver,
										const Point3d& normal,
										const Point3d& point) const
{
	// расчёт коор-т в СК наблюдателя
	const Point3d p_pr = point - normal*Point3d::DotProduct(normal, point);
	return Point3d(Point3d::DotProduct(p_pr, hor),
				   Point3d::DotProduct(p_pr, ver), 0);
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
		ver = Point3d::CrossProduct(normal, hor);
	}
}

void Handler::ComputeLengthIndices(const Beam &beam, BeamInfo &info) const
{
	auto *lens = info.opticalLengths;

	Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.vertices[0]) - info.projectedCenter;
	Point3d p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.vertices[1]) - info.projectedCenter;
	Point3d p3 = ChangeCoordinateSystem(info.horAxis, info.verAxis, info.normald,
										beam.vertices[2]) - info.projectedCenter;

	double den = p1.x*(p2.y - p3.y) - p2.x*(p1.y + p3.y) + p3.x*(p1.y - p2.y);

	info.lenIndices.z = (lens[0]*(p2.x*p3.y - p3.x*p2.y) -
			lens[1]*(p1.x*p3.y + p3.x*p1.y) +
			lens[2]*(p1.x*p2.y - p2.x*p1.y)) / den;

	info.lenIndices.x = (lens[0]*(p2.y - p3.y) -
			lens[1]*(p1.y + p3.y) +
			lens[2]*(p1.y - p2.y)) / den;

	info.lenIndices.y = -(lens[0]*(p2.x - p3.x) -
			lens[1]*(p1.x + p3.x) +
			lens[2]*(p1.x - p2.x)) / den;
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
	Point3d k_k0 = -direction +
			Point3d(dir.coordinates[0], dir.coordinates[1], dir.coordinates[2]);

	Point3d	pt_proj = ChangeCoordinateSystem(info.horAxis, info.verAxis,
											 info.normald, k_k0);

	const complex A(pt_proj.x, info.lenIndices.x*m_riIm);
	const complex B(pt_proj.y, info.lenIndices.y*m_riIm);

	if (abs(A) < m_eps2 && abs(B) < m_eps2)
	{
		return m_invComplWave * info.area;
	}

	complex s(0, 0);

	int begin, startIndex, endIndex;

	if (info.order)
	{
		begin = 0;
		startIndex = beam.nVertices-1;
		endIndex = -1;
	}
	else
	{
		begin = beam.nVertices-1;
		startIndex = 0;
		endIndex = beam.nVertices;
	}

	Point3d p1 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.vertices[begin]) - info.projectedCenter;
	Point3d p2;

	if (abs(B) > abs(A))
	{
		for (int i = startIndex; i != endIndex; i = info.order ? i-1 : i+1)
		{
			p2 = ChangeCoordinateSystem(info.horAxis, info.verAxis,
										info.normald, beam.vertices[i]) - info.projectedCenter;

			if (fabs(p1.x - p2.x) < m_eps1)
			{
				p1 = p2;
				continue;
			}

			const double ai = (p1.y - p2.y)/(p1.x - p2.x);
			complex Ci = A + ai*B;

			complex tmp;

			if (abs(Ci) < m_eps1)
			{
				double mul = p2.x*p2.x - p1.x*p1.x;
				tmp = complex(-m_wi2*real(Ci)*mul/2.0,
							  m_wavenumber*(p2.x - p1.x) + m_wi2*imag(Ci)*mul/2.0);
			}
			else
			{
				double kReCi = m_wavenumber*real(Ci);
				double kImCi = -m_wavenumber*imag(Ci);
				tmp = (exp_im(kReCi*p2.x)*exp(kImCi*p2.x) -
					   exp_im(kReCi*p1.x)*exp(kImCi*p1.x))/Ci;
			}

			const double bi = p1.y - ai*p1.x;
			double kBi = m_wavenumber*bi;
			s += exp_im(kBi*real(B)) * tmp * exp(-kBi*imag(B));

			p1 = p2;
		}

		s /= B;
	}
	else
	{
		for (int i = startIndex; i != endIndex; i = info.order ? i-1 : i+1)
		{
			p2 = ChangeCoordinateSystem
					(info.horAxis, info.verAxis, info.normald, beam.vertices[i])
					- info.projectedCenter;

			if (fabs(p1.y - p2.y) < m_eps1)
			{
				p1 = p2;
				continue;
			}

			const double ci = (p1.x - p2.x)/(p1.y - p2.y);
			const complex Ei = A*ci + B;

			complex tmp;

			if (abs(Ei) < m_eps1)
			{
				double mul = p2.y*p2.y - p1.y*p1.y;
				tmp = complex(-m_wi2*real(Ei)*mul/2.0, // REF вынести m_wi2*mul/2.0
							  m_wavenumber*(p2.y - p1.y) + m_wi2*imag(Ei)*mul/2.0);
			}
			else
			{
				double kReEi = m_wavenumber*real(Ei);
				double kImEi = -m_wavenumber*imag(Ei);
				tmp = (exp_im(kReEi*p2.y)*exp(kImEi*p2.y) -
					   exp_im(kReEi*p1.y)*exp(kImEi*p1.y))/Ei;
			}

			const double di = p1.x - ci*p1.y;
			double kDi = m_wavenumber*di;
			s += exp_im(kDi*real(A)) * exp(-kDi*imag(A)) * tmp;

			p1 = p2;
		}

		s /= -A;
	}

	return m_complWave * s * exp(m_absMag*info.lenIndices.z); // REF вынести
}

complex Handler::DiffractIncline(const BeamInfo &info, const Beam &beam,
								 const Point3d &direction) const
{
	const Point3f &dir = beam.direction;
	Point3d k_k0 = -direction + Point3d(dir);

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

	int begin, startIndex, endIndex;

	if (info.order)
	{
		begin = 0;
		startIndex = beam.nVertices-1;
		endIndex = -1;
	}
	else
	{
		begin = beam.nVertices-1;
		startIndex = 0;
		endIndex = beam.nVertices;
	}

	Point3d p1 = ChangeCoordinateSystem
			(info.horAxis, info.verAxis, info.normald, beam.vertices[begin]) -
			info.projectedCenter;
	Point3d p2;

	if (absB > absA)
	{
		for (int i = startIndex; i != endIndex; i = info.order ? i-1 : i+1)
		{
			p2 = ChangeCoordinateSystem
					(info.horAxis, info.verAxis, info.normald, beam.vertices[i])
					- info.projectedCenter;

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
		for (int i = startIndex; i != endIndex; i = info.order ? i-1 : i+1)
		{
			p2 = ChangeCoordinateSystem
					(info.horAxis, info.verAxis, info.normald, beam.vertices[i])
					- info.projectedCenter;

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
