#include "Distributions.h"

template<class Distribution1, class Distribution2>
class Mixture : public IDistribution, public IPersistent
{
private:

	Distribution1* D1;

	Distribution2* D2;

	double p;

public:

	Mixture(Distribution1* D1, Distribution2* D2, double p) : D1(D1), D2(D2), p(p) {};

	Distribution1* get_component1();

	Distribution2* get_component2();

	double get_p() const;

	void set_p(const double p);

	void load_from_file(ifstream& file) override;

	void save_to_file(ofstream& file) override;

	double f(const double x) const override;

	double expected_value() const override;

	double variance() const override;

	double asymmetry() const override;

	double kurtosis() const override;

	double random_var() const override;

	vector<double> generate_sequence(const int n) const override;

	vector<pair<double, double>> generate_table_of_values(const int n, const vector<double>& x_s = {}) const override;
};

template <class Distribution1, class Distribution2>
Distribution1* Mixture<Distribution1, Distribution2>::get_component1()
{
	return D1;
}

template <class Distribution1, class Distribution2>
Distribution2* Mixture<Distribution1, Distribution2>::get_component2()
{
	return D2;
}


template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::get_p() const
{
	return p;
}

template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::set_p(const double p)
{
	if (p < 0 || p > 1)
		throw invalid_argument("Значение параметра некорректно");
	this->p = p;
}

template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::load_from_file(ifstream& file)
{
	string filename;
	//file.open("mixture_params.txt");
	ifstream file1;
	ifstream file2;

	cout << "Введите имя файла с параметром смеси распределений: ";
	cin >> filename;

	file.open(filename);
	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	file >> p;
	get_component1()->load_from_file(file1);
	get_component2()->load_from_file(file2);

	file.close();
}

template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::save_to_file(ofstream& file)
{
	ofstream file1;
	ofstream file2;
	file.open("mixture_params.txt");

	file << p;
	get_component1()->save_to_file(file1);
	get_component2()->save_to_file(file2);

	file.close();

	cout << "Параметр смеси сохранен в файл mixture_params.txt" << endl;

}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::f(const double x) const
{
	return (1 - p) * D1->f(x) + p * D2->f(x);
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::expected_value() const
{
	return (1 - p) * D1->expected_value() + p * D2->expected_value();
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::variance() const
{
	return (1 - p) * (pow(D1->expected_value(), 2) + D1->variance()) +
		p * (pow(D2->expected_value(), 2) + D2->variance()) -
		pow(expected_value(), 2);
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::asymmetry() const
{
	return ((1 - p) * (pow((D1->expected_value() - expected_value()), 3) + 3 * (D1->expected_value() - expected_value()) * D1->variance() + pow(D1->variance(), 3 / 2) * D1->asymmetry()) +
		p * (pow((D2->expected_value() - expected_value()), 3) + 3 * (D2->expected_value() - expected_value()) * D2->variance() + pow(D2->variance(), 3 / 2) * D2->asymmetry())) /
		pow(variance(), 3 / 2);
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::kurtosis() const
{
	return ((1 - p) * (pow((D1->expected_value() - expected_value()), 4) + 6 * D1->variance() * pow((D1->expected_value() - expected_value()), 2) +
		4 * (D1->expected_value() - expected_value()) * pow(D1->variance(), 3 / 2) * D1->asymmetry() + pow(D1->variance(), 2) * D1->kurtosis()) +
		p * (pow((D2->expected_value() - expected_value()), 4) + 6 * D2->variance() * pow((D2->expected_value() - expected_value()), 2) +
			4 * (D2->expected_value() - expected_value()) * pow(D2->variance(), 3 / 2) * D2->asymmetry() + pow(D2->variance(), 2) * D2->kurtosis()) - 3) /
		pow(variance(), 2);
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::random_var() const
{
	random_device rd;
	default_random_engine gen(rd());
	uniform_real_distribution<> d(0, 1);

	double r = d(gen);

	if (r > p)
		return D1->random_var();
	else
		return D2->random_var();
}

template <class Distribution1, class Distribution2>
vector<double> Mixture<Distribution1, Distribution2>::generate_sequence(const int n) const
{
	vector<double> sequence;
	for (int i = 0; i < n; i++)
	{
		double x = random_var();
		sequence.push_back(x);
	}

	sort(sequence.begin(), sequence.end());
	return sequence;
}

template <class Distribution1, class Distribution2>
vector<pair<double, double>> Mixture<Distribution1, Distribution2>::generate_table_of_values(const int n, const vector<double>& x_s) const
{
	vector<double> sequence;
	vector<pair<double, double>> table;

	if (x_s.empty())
		sequence = generate_sequence(n);
	else
		sequence = x_s;

	for (const double& x : sequence)
	{
		double y = f(x);
		table.push_back(make_pair(x, y));
	}

	return table;
}
