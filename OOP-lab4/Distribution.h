#pragma once

__interface IDistribution
{
public:
	double virtual density(const double x) const = 0;

	double virtual M_Ksi() const = 0;

	double virtual D_Ksi() const = 0;

	double virtual asymmetry() const = 0;

	double virtual kurtosis() const = 0;

	double virtual algorithm() const = 0;

	vector<double> virtual selection(const int n) const = 0;

	vector<pair<double, double>> virtual generate_pair(const int n, const vector<double>& x_selection = {}) const = 0;
};

__interface IPersistent
{
public:
	void virtual load_file(ifstream& file) = 0;

	void virtual save_to_file(ofstream& file) = 0;
};
