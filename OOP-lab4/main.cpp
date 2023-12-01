#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "function.h"
#include "Distribution.h"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <string>
using namespace std;

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");
	int n = 1000;
	int result = Catch::Session().run(argc, argv);
	return result;
}