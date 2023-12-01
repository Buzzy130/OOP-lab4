#pragma once
#include "Lib.h"

__interface IDistribution
{
public:
	double virtual f(const double x) const = 0;

	double virtual expected_value() const = 0;

	double virtual variance() const = 0;

	double virtual asymmetry() const = 0;

	double virtual kurtosis() const = 0;

	double virtual random_var() const = 0;

	vector<double> virtual generate_sequence(const int n) const = 0;

	vector<pair<double, double>> virtual generate_table_of_values(const int n, const vector<double>& x_s = {}) const = 0;
};

__interface IPersistent
{
public:
	void virtual load_from_file(ifstream& file) = 0;

	void virtual save_to_file(ofstream& file) = 0;
};