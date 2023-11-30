#pragma once

__interface IDistribution
{
public:
	double virtual Huber_D(const double x) const = 0;

	double virtual M_Ksi() const = 0;

	double virtual D_Ksi() const = 0;

	double virtual asymmetry() const = 0;

	double virtual kurtosis() const = 0;

	double virtual random_var() const = 0;

	vector<double> virtual generate_sequence(const int n) const = 0;

	vector<pair<double, double>> virtual generate_table_of_values(const int n, const vector<double>& x_s = {}) const = 0;
};

__interface IPersistent
{
public:
	void virtual load_file(ifstream& file) = 0;

	void virtual save_to_file(ofstream& file) = 0;
};
