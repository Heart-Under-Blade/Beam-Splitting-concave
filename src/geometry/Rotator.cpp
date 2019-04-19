#include "Rotator.h"

void LocalRotator::SetRotationAngle(const Orientation &angle)
{
	double cosB, cosG,
			sinB, sinG;

	sincos(angle.zenith,  &sinB, &cosB);
	sincos(angle.azimuth, &sinG, &cosG);

	m_rotationMatrix[0][0] = cosB*cosG - 0;
	m_rotationMatrix[1][0] = sinG;
	m_rotationMatrix[2][0] = -sinB*cosG;

	m_rotationMatrix[0][1] = -(cosB*sinG + 0);
	m_rotationMatrix[1][1] = cosG;
	m_rotationMatrix[2][1] = sinB*sinG;

	m_rotationMatrix[0][2] = sinB;
	m_rotationMatrix[1][2] = 0;
	m_rotationMatrix[2][2] = cosB;
}

void GlobalRotator::SetRotationAngle(const Orientation &angle)
{
	double cosF, cosT, cosP,
			sinF, sinT, sinP;

	sincos(0, &sinP, &cosP);
	sincos(angle.zenith, &sinF, &cosF);
	sincos(angle.azimuth, &sinT, &cosT);

	double sinFsinT = sinF*sinT;
	double cosFsinT = cosF*sinT;

	m_rotationMatrix[0][0] = cosT*cosP;
	m_rotationMatrix[1][0] = sinFsinT*cosP + cosF*sinP;
	m_rotationMatrix[2][0] = -cosFsinT*cosP + sinF*sinP;

	m_rotationMatrix[0][1] = -cosT*sinP;
	m_rotationMatrix[1][1] = -sinFsinT*sinP + cosF*cosP;
	m_rotationMatrix[2][1] = cosFsinT*sinP + sinF*cosP;

	m_rotationMatrix[0][2] = sinT;
	m_rotationMatrix[1][2] = -sinF*cosT;
	m_rotationMatrix[2][2] = cosF*cosT;
}

void Rotator::RotatePoint(const Point3f &point, Point3f &result)
{
	result.coordinates[0] = point.coordinates[0] * m_rotationMatrix[0][0] +
			point.coordinates[1] * m_rotationMatrix[0][1] +
			point.coordinates[2] * m_rotationMatrix[0][2];

	result.coordinates[1] = point.coordinates[0] * m_rotationMatrix[1][0] +
			point.coordinates[1] * m_rotationMatrix[1][1] +
			point.coordinates[2] * m_rotationMatrix[1][2];

	result.coordinates[2] = point.coordinates[0] * m_rotationMatrix[2][0] +
			point.coordinates[1] * m_rotationMatrix[2][1] +
			point.coordinates[2] * m_rotationMatrix[2][2];
}
