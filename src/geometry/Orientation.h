#pragma once

#include <cmath>

/**
 * @brief Spatial orientation with two angled coordinates
 */
class Orientation
{
public:
    double zenith;
    double azimuth;

    Orientation();
    Orientation(double b, double g);

    /**
     * @brief Convert angles to radians
     */
    void ToRadians();

    Orientation ToRadians() const;

    /**
     * @brief Convert angles to degrees
     */
    void ToDegrees();

    static double DegToRad(double deg);
    static double RadToDeg(double rad);
};

class OrientationRange
{
public:
	int nZenith;
	int nAzimuth;
	Orientation from;
	Orientation to;
	Orientation step;
};
