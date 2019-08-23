#pragma once

#include "Column.h"

class RegularColumn : public Column
{
public:
	RegularColumn(const complex &refrIndex, const Size &size);
};
