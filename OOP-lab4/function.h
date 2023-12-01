#include "Distribution.h"

class Primary : public IDistribution, public IPersistent
{
private:

	double v;

	double k;

	double scale;

	double shift;

public:

	Primary(double v = 1, double scale = 1, double shift = 0);

	Primary(ifstream& file);

	void save_to_file(ofstream& file) override;

	void load_from_file(ifstream& file) override;

	double get_v() const;

	void set_v(const double v);

	double get_k() const;

	double get_scale() const;

	void set_scale(const double scale);

	double get_shift() const;

	void set_shift(const double shift);

	double f(const double x) const override;

	double phi(const double x) const;

	double phi_lower(const double x) const;

	double expected_value() const override;

	double variance() const override;

	double asymmetry() const override;

	double kurtosis() const override;

	double P() const;

	double K(const double v) const;

	double random_var() const override;

	vector<double> generate_sequence(const int n) const override;

	vector<pair<double, double>> generate_table_of_values(const int n, const vector<double>& x_s = {}) const override;
};
//-----------------------------------------emprical-------------------------------------------------//
class Empirical : public IDistribution, public IPersistent
{
private:

	vector<double> x_s;

	vector<double> f_s;

	int n = 0;

	int k = 0;

public:

	Empirical(const IDistribution* D, int _n, int _k);

	Empirical(const Empirical* EM);

	Empirical(const int _n, const int _k);

	Empirical(const vector<double>& x_s);

	Empirical(ifstream& file);

	~Empirical();

	Empirical& operator=(const Empirical& EM);

	double random_var() const override;

	vector<double> generate_sequence(const int n) const override;

	vector<double> generate_values() const;

	vector<pair<double, double>> generate_table_of_values(const int n, const vector<double>& x_s = {}) const override;

	vector<double> get_x_s() const;

	vector<double> get_f_s() const;

	int get_n() const;

	int get_k() const;

	void save_to_file(ofstream& file) override;

	void load_from_file(ifstream& file) override;

	double f(const double x) const override;

	double expected_value() const override;

	double variance() const override;

	double asymmetry() const override;

	double kurtosis() const override;
};