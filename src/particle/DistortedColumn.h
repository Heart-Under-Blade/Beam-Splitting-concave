#pragma once

#include "Column.h"

/**
 * @brief The Hexagon class
 * The prism particle with 6 number of side facets.
 */
class DistortedColumn : public Column
{
public:
	DistortedColumn();
	DistortedColumn(const complex &refrIndex, double diameter, double height,
					   double angle);

private:
	void DistortBases(double angle);
};

