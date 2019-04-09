#include "Orientation.h"

Orientation::Orientation()
{
}

Orientation::Orientation(double b, double g)
{
    zenith = b;
    azimuth = g;
}

void Orientation::ToRadians()
{
    zenith = DegToRad(zenith);
    azimuth = DegToRad(azimuth);
}

Orientation Orientation::ToRadians() const
{
    Orientation angle = *this;
    angle.ToRadians();
    return angle;
}

void Orientation::ToDegrees()
{
    zenith = RadToDeg(zenith);
    azimuth = RadToDeg(azimuth);
}

double Orientation::DegToRad(double deg)
{
    return (deg*M_PI)/180;
}

double Orientation::RadToDeg(double rad)
{
    return (rad*180)/M_PI;
}
