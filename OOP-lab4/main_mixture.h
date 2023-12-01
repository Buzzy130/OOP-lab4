#pragma once
#include "mixture.cpp"
#include "function.h"

template<class Distribution1, class Distribution2>
void final_mixture(Mixture<Distribution1, Distribution2>* MX);

void next_mixture();