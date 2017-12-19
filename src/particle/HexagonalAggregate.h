#pragma once

#include "Hexagonal.h"

class HexagonalAggregate : public Hexagonal
{
public:
	HexagonalAggregate(const complex &refrIndex, double diameter, double height,
					   int particleNumber);

	bool IsComplicated() const override;

protected:
	void SetFacetParams() override;

protected:
//	static const double PART_OFFSET = 1;
};
