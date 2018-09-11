#include "Particle.h"
#include <fstream>
#include <iostream>
#include <cstring>
#include "global.h"

Particle::Particle()
{
	SetFacetIndices();
}

Particle::Particle(int nFacets, const complex &refrIndex, bool isNonConvex)
	: m_refractiveIndex(refrIndex),
	  m_isNonConvex(isNonConvex)
{
	nElems = nFacets;
	SetFacetIndices();
}

void Particle::SetFacetIndices()
{
	for (int i = 0; i < MAX_FACET_NUM; ++i)
	{
		elems[i].actual.index = i;
	}
}

Facet *Particle::GetActualFacet(int i)
{
	return &elems[i].actual;
}

void Particle::SetFromFile(const std::string &filename)
{
	std::ifstream pfile(filename, std::ios::in);

	if (!pfile.is_open())
	{
		std::cerr << "File \"" << filename << "\" is not found" << std::endl;
		throw std::exception();
	}

	const int bufSize = 1024;
	char *buff = (char*)malloc(sizeof(char) * bufSize);

	nElems = 0;
	Facet &facet = elems[nElems++].origin;

	char *ptr, *trash;

	pfile.getline(buff, bufSize);
	m_isNonConvex = strtol(buff, &trash, 10);

	pfile.getline(buff, bufSize);
	isAggregated = strtol(buff, &trash, 10);

	// read symmetry params
	{
		pfile.getline(buff, bufSize);

		ptr = strtok(buff, " ");
		m_symmetry.beta = Angle::DegToRad(strtod(ptr, &trash));

		ptr = strtok(NULL, " ");
		m_symmetry.gamma = Angle::DegToRad(strtod(ptr, &trash));
	}

	pfile.getline(buff, bufSize); // skip empty line

	while (!pfile.eof()) // read vertices of facets
	{
		pfile.getline(buff, bufSize);
		ptr = strtok(buff, " ");

		if (strlen(buff) == 0)
		{
			facet = elems[nElems++].origin;
			continue;
		}

		int c_i = 0;

		while (ptr != NULL)
		{
			facet.arr[facet.nVertices].point[c_i++] = strtod(ptr, &trash);
			ptr = strtok(NULL, " ");
		}

		++(facet.nVertices);
	}

	pfile.close();

	// correction of number of facet
	if (elems[nElems-1].origin.nVertices == 0)
	{
		--nElems;
	}

	SetDefaultNormals();
	Reset();
	SetDefaultCenters();

	if (IsNonConvex() || isAggregated)
	{
		for (int i = 0; i < nElems; ++i)
		{
			elems[i].origin.isOverlayedIn = false;
			elems[i].origin.isOverlayedOut = false;
			elems[i].actual.isOverlayedIn = false;
			elems[i].actual.isOverlayedOut = false;
		}
	}
}

void Particle::Rotate(const Angle &angle)
{
	rotAngle = angle;
	SetRotateMatrix();

	// REF: слить всё в один цикл
	for (int i = 0; i < nElems; ++i)
	{
		auto &facet = elems[i];

		for (int j = 0; j < facet.origin.nVertices; ++j)
		{
			RotatePoint(facet.origin.arr[j], facet.actual.arr[j]);
		}

		// centers
		RotatePoint(facet.origin.center, facet.actual.center);
	}

	RotateNormals();
}

void Particle::Fix()
{
	for (int i = 0; i < nElems; ++i)
	{
		for (int j = 0; j < elems[i].actual.nVertices; ++j)
		{
			elems[i].origin.arr[j] = elems[i].actual.arr[j];
		}
	}
}

void Particle::Concate(const std::vector<Particle> &parts)
{
	int i = 0;
	nElems = 0;

	for (const Particle &part : parts)
	{
		nElems += part.nElems;

		for (int j = 0; j < part.nElems; ++j)
		{
			elems[i++].origin = part.elems[j].actual;
		}
	}

	isAggregated = true;
}

double Particle::GetRotationRadius() const
{
	Point3f p0(0, 0, 0);

	double radius = 0;

	for (int i = 0; i < nElems; ++i)
	{
		const Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			Point3f v_len = facet.arr[j] - p0;
			double len = Length(v_len);

			if (len > radius)
			{
				radius = len;
			}
		}
	}

	return radius;
}

const complex &Particle::GetRefractiveIndex() const
{
	return m_refractiveIndex;
}

const Symmetry &Particle::GetSymmetry() const
{
	return m_symmetry;
}

void Particle::Move(float dx, float dy, float dz)
{
	for (int i = 0; i < nElems; ++i)
	{
		for (int j = 0; j < elems[i].origin.nVertices; ++j)
		{
			elems[i].actual.arr[j] = elems[i].origin.arr[j] + Point3f(dx, dy, dz);
		}
	}
}

bool Particle::IsNonConvex() const
{
	return m_isNonConvex;
}

void Particle::Output()
{
	std::ofstream M("particle.dat", std::ios::out);

	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;

		for (int j = 0; j < facet.nVertices; ++j)
		{
			Point3f p = facet.arr[j];
			M << p.point[0] << ' '
							<< p.point[1] << ' '
							<< p.point[2] << ' '
							<< i ;
			M << std::endl;
		}

		M << std::endl << std::endl;;
	}

	M.close();
}

void Particle::SetRefractiveIndex(const complex &value)
{
	m_refractiveIndex = value;
}

void Particle::SetDefaultNormals()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].origin.SetNormal();
	}
}

void Particle::Reset()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].actual = elems[i].origin;
	}
}

void Particle::SetDefaultCenters()
{
	for (int i = 0; i < nElems; ++i)
	{
		elems[i].origin.SetCenter();
	}
}

void Particle::SetRotateMatrix()
{
	double cosA, cosB, cosG,
			sinA, sinB, sinG;

	sincos(rotAngle.alpha, &sinA, &cosA);
	sincos(rotAngle.beta,  &sinB, &cosB);
	sincos(rotAngle.gamma, &sinG, &cosG);

	double cosAcosB = cosA*cosB;
	double sinAcosG = sinA*cosG;
	double sinAsinG = sinA*sinG;

	m_rotMatrix[0][0] = cosAcosB*cosG - sinAsinG;
	m_rotMatrix[1][0] = sinAcosG*cosB + cosA*sinG;
	m_rotMatrix[2][0] = -sinB*cosG;

	m_rotMatrix[0][1] = -(cosAcosB*sinG + sinAcosG);
	m_rotMatrix[1][1] = cosA*cosG - sinAsinG*cosB;
	m_rotMatrix[2][1] = sinB*sinG;

	m_rotMatrix[0][2] = cosA*sinB;
	m_rotMatrix[1][2] = sinA*sinB;
	m_rotMatrix[2][2] = cosB;
}

void Particle::RotateNormals()
{
	for (int i = 0; i < nElems; ++i)
	{
		auto &facet = elems[i];
		RotatePoint(facet.origin.in_normal, facet.actual.in_normal);
	}

	SetDParams();

	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;
		facet.ex_normal = -facet.in_normal;
		facet.ex_normal.d_param = -facet.in_normal.d_param;
	}
}

void Particle::SetDParams()
{
	for (int i = 0; i < nElems; ++i)
	{
		Facet &facet = elems[i].actual;
		double d = DotProduct(facet.arr[0], facet.in_normal);
		facet.in_normal.d_param = -d;
	}
}

void Particle::RotatePoint(const Point3f &point, Point3f &result)
{
	result.cx = point.cx*m_rotMatrix[0][0] + point.cy*m_rotMatrix[0][1] + point.cz*m_rotMatrix[0][2];
	result.cy = point.cx*m_rotMatrix[1][0] + point.cy*m_rotMatrix[1][1] + point.cz*m_rotMatrix[1][2];
	result.cz = point.cx*m_rotMatrix[2][0] + point.cy*m_rotMatrix[2][1] + point.cz*m_rotMatrix[2][2];
}

void Particle::SetSymmetry(double beta, double gamma, double alpha)
{
	m_symmetry.beta = beta;
	m_symmetry.gamma = gamma;
	m_symmetry.alpha = alpha;
}
