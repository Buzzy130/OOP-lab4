#pragma once
#include "Distributions.h"

class Huber : public IDistribution, public IPersistent
{
private:

	double v;

	double k;

	double scale;

	double shift;

public:

	Huber(double v = 1, double scale = 1, double shift = 0);

	Huber(ifstream& file);

	void save_to_file(ofstream& file) override;

	void load_from_file(ifstream& file) override;

	double get_v() const;

	void set_v(const double v);

	double get_k() const;

	double get_scale() const;

	void set_scale(const double scale);

	double get_shift() const;

	void set_shift(const double shift);

	double density(const double x) const override;

	double phi(const double x) const;

	double phi_lower(const double x) const;

	double M_Ksi() const override;

	double D_Ksi() const override;

	double asymmetry() const override;

	double kurtosis() const override;

	double P() const;

	double K(const double v) const;

	double algorithm() const override;

	vector<double> selection(const int n) const override;

	vector<pair<double, double>> generate_pair(const int n, const vector<double>& x_selection = {}) const override;
};
//========================================================================//

class Empirical : public IDistribution, public IPersistent
{
private:

	vector<double> x_selection;

	vector<double> f_selection;

	int n = 0;

	int k = 0;

public:

	Empirical(const IDistribution* D, int _n, int _k);

	Empirical(const Empirical* EM);

	Empirical(const int _n, const int _k);

	Empirical(const vector<double>& x_selection);

	Empirical(ifstream& file);

	~Empirical();

	Empirical& operator=(const Empirical& EM);

	double algorithm() const override;

	vector<double> selection(const int n) const override;

	vector<double> generate_f_selection() const;

	vector<pair<double, double>> generate_pair(const int n, const vector<double>& x_selection = {}) const override;

	vector<double> get_x_selection() const;

	vector<double> get_f_selection() const;

	int get_n() const;

	int get_k() const;

	void save_to_file(ofstream& file) override;

	void load_from_file(ifstream& file) override;

	double density(const double x) const override;

	double M_Ksi() const override;

	double D_Ksi() const override;

	double asymmetry() const override;

	double kurtosis() const override;
};