#pragma once

#include "Column.h"

class HexagonalAggregate : public Column
{
public:
	HexagonalAggregate(const Size &size, int particleNumber);
protected:
	void SetFacetParams() override;

protected:
//	static const double PART_OFFSET = 1;
};
