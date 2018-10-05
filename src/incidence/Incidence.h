#pragma once

class Beam;
class Splitting;

class Incidence
{
public:
	virtual void ComputeDirections(const Beam &beam,
									Splitting &splitter) const {};

    virtual void ComputeJonesMatrices(const Beam &beam,
                                      Splitting &splitter) const {};

    void ComputeOpticalPaths(const Beam &beam, Splitting &splitter) const;
};
