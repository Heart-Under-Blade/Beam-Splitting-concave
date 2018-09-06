#pragma once

#include "Column.h"

class HexagonalAggregate : public Column
{
public:
	HexagonalAggregate(const complex &refrIndex, const Size &size,
					   int particleNumber);
protected:
	void SetFacetParams() override;

protected:
//	static const double PART_OFFSET = 1;
};
