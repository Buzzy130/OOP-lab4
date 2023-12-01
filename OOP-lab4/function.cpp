#include "function.h"
#include "Distribution.h"

Primary::Primary(double v, double scale, double shift)
{
	if (v <= 0 || scale <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров не корректны");

	this->v = v;
	this->k = K(v);
	this->scale = scale;
	this->shift = shift;
}

Primary::Primary(ifstream& file)
{
	load_from_file(file);
}

void Primary::save_to_file(ofstream& file)
{
	string filename;
	cout << "Введите имя файла, куда следует записать параметры основного распределения: ";
	cin >> filename;

	file.open(filename);
	file << this->v << endl << this->scale << endl << this->shift;
	file.close();

	cout << "Параметры основного распределения сохранены в файл " + filename << endl;
}

void Primary::load_from_file(ifstream& file)
{
	string filename;
	double v, scale, shift;

	cout << "Введите имя файла, откуда следует считать параметры основного распределения: ";
	cin >> filename;

	file.open(filename);

	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	file >> v >> scale >> shift;

	file.close();

	if (v <= 0 || scale <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров некорректны");

	this->v = v;
	this->scale = scale;
	this->shift = shift;
}


double Primary::get_v() const
{
	return this->v;
}

double Primary::get_k() const
{
	return this->k;
}


void Primary::set_v(double v)
{
	if (v <= 0)
		throw invalid_argument("Ошибка: v <= 0");

	this->v = v;
	this->k = K(this->v);
}

double Primary::get_scale() const
{
	return this->scale;
}

void Primary::set_scale(const double scale)
{
	if (scale <= 0)
		throw invalid_argument("Ошибка: scale <= 0");
	this->scale = scale;
}

double Primary::get_shift() const
{
	return this->shift;
}

void Primary::set_shift(const double shift)
{
	this->shift = shift;
}

double Primary::f(const double x) const
{
	return (1. / (sqrt(2. * M_PI) * this->k) * (abs((x - this->shift) / this->scale) <= this->v ? exp(-pow((x - this->shift) / this->scale, 2.) / 2.) : exp(pow(this->v, 2.) / 2. - this->v * abs((x - this->shift) / this->scale)))) / this->scale;
}

double Primary::phi(double x) const
{
	return 0.5 * (1. + erf(x / sqrt(2.)));
}

double Primary::phi_lower(double x) const
{
	return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}

double Primary::expected_value() const
{
	return this->shift;
}

double Primary::variance() const
{
	return 1. + 2. * phi_lower(this->v) * (pow(this->v, 2.) + 2.) / (pow(this->v, 3.) * this->k);
}

double Primary::asymmetry() const
{
	return 0.;
}

double Primary::kurtosis() const
{
	return (3. * (2. * phi(this->v) - 1.) + 2. * phi_lower(this->v) * (24. / pow(this->v, 5.) + 24. / pow(this->v, 3.) + 12. / this->v + this->v)) / (pow(variance(), 2.) * this->k) - 3.;
}

double Primary::P() const
{
	return (2. * phi(this->v) - 1.) / this->k;
}

double Primary::K(const double v) const
{
	return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}

double Primary::random_var() const
{
	random_device rd;
	default_random_engine gen(rd());
	uniform_real_distribution<> d(0, 1);

	//шаг 1
	double r1 = d(gen);


	if (r1 <= P())
	{
		//шаг 2
		double r2, r3, x1;

		do {
			r2 = d(gen);
			r3 = d(gen);
			x1 = sqrt(-2 * log(r2)) * cos(2 * M_PI * r3);
			//double x1 = sqrt(-2 * log(r2)) * sin(2 * M_PI * r3)
		} while (!(-this->v <= x1 && x1 <= this->v)); //шаг 3

		return x1 * this->scale + this->shift;
	}
	else
	{
		//шаг 4
		double r4 = d(gen);
		double x2 = this->v - log(r4) / this->v;

		//шаг 5
		return r1 < (1 + P()) / 2 ? x2 * this->scale + this->shift : -x2 * this->scale + this->shift;
	}
}

vector<double> Primary::generate_sequence(const int n) const
{
	vector<double> res;

	for (int i = 0; i < n; i++)
	{
		double x = random_var();
		res.push_back(x);
	}

	sort(res.begin(), res.end());

	return res;
}

vector<pair<double, double>> Primary::generate_table_of_values(const int n, const vector<double>& x_s) const
{
	vector<pair<double, double>> res;
	vector<double> sequence;

	if (x_s.empty())
		sequence = generate_sequence(n);
	else
		sequence = x_s;

	for (const double& x : sequence)
	{
		double y = f(x);
		res.push_back(make_pair(x, y));
	}

	return res;
}


//---------------------------------------mixture--------------------------------------------------//


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
		throw invalid_argument("Çíà÷åíèå ïàðàìåòðà íåêîððåêòíî");
	this->p = p;
}

template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::load_from_file(ifstream& file)
{
	string filename;
	//file.open("mixture_params.txt");
	ifstream file1;
	ifstream file2;

	cout << "Ââåäèòå èìÿ ôàéëà ñ ïàðàìåòðîì ñìåñè ðàñïðåäåëåíèé: ";
	cin >> filename;

	file.open(filename);
	if (!file)
		throw runtime_error("Îøèáêà: íå óäàëîñü îòêðûòü ôàéë");

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

	cout << "Ïàðàìåòð ñìåñè ñîõðàíåí â ôàéë mixture_params.txt" << endl;

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

//------------------------------------------------Emprical--------------------------------------------------//
Empirical::Empirical(const IDistribution* D, int _n, int _k) :
	n(_n > 1 ? _n : throw invalid_argument("Некорректный аргумент")), k(_k > 2 ? _k : (int)floor(log2(_n)) + 1)
{
	x_s = D->generate_sequence(_n);
	f_s = generate_values();
}

Empirical::Empirical(const Empirical* EM) : n(EM->n > 1 ? EM->n : throw invalid_argument("Некорректный аргумент")), k(EM->k > 2 ? EM->k : (int)floor(log2(EM->n) + 1)), x_s(EM->x_s), f_s(EM->f_s)
{
	f_s = generate_values();
}

Empirical::Empirical(const int _n, const int _k) : n(_n > 1 ? _n : throw invalid_argument("Некорректный аргумент")), k(_k > 2 ? _k : (int)floor(log2(_n) + 1))
{
	x_s = generate_sequence(n);
	f_s = generate_values();
}

Empirical::Empirical(const vector<double>& x_s) : n(x_s.size()), k((int)floor(log2(x_s.size())) + 1)
{
	this->x_s = x_s;
	f_s = generate_values();
}

Empirical::Empirical(ifstream& file)
{
	load_from_file(file);
}

Empirical::~Empirical()
{
	x_s.clear();
	f_s.clear();
}

Empirical& Empirical::operator=(const Empirical& EM)
{
	if (this == &EM)
		return *this;

	x_s = EM.x_s;
	f_s = EM.f_s;
	n = EM.n;
	k = EM.k;
	return *this;
}

double Empirical::random_var() const
{
	vector<double> intervals;
	vector<double> densities;

	int k = get_k();
	double delta = (1.0 - 0.0) / (double)k;

	double point = 0.0;
	while (point < 1.0)
	{
		intervals.push_back(point);
		point += delta;
	}
	//intervals.push_back(1.0);

	for (int i = 0; i < k; i++)
		densities.push_back(delta);

	random_device rd;
	default_random_engine gen(rd());
	piecewise_constant_distribution<> d(intervals.begin(), intervals.end() - 1, densities.begin());

	return d(gen);

}

vector<double> Empirical::generate_sequence(const int n) const
{
	vector<double> result;

	for (int i = 0; i < n; i++)
		result.push_back(random_var());

	sort(result.begin(), result.end());

	return result;
}


vector<double> Empirical::generate_values() const
{
	vector<double> result;

	for (const double& x : x_s)
		result.push_back(f(x));

	return result;
}

vector<pair<double, double>> Empirical::generate_table_of_values(const int n, const vector<double>& x_s) const
{
	vector<pair<double, double >> result;

	for (int i = 0; i < n; i++)
		result.push_back(make_pair(x_s[i], f_s[i]));

	return result;
}

vector<double> Empirical::get_x_s() const
{
	return x_s;
}

vector<double> Empirical::get_f_s() const
{
	return f_s;
}

int Empirical::get_n() const
{
	return n;
}

int Empirical::get_k() const
{
	return k;
}

void Empirical::save_to_file(ofstream& file)
{
	vector<pair<double, double>> pairs = generate_table_of_values(n);

	ofstream file_sequence;
	ofstream file_values;
	ofstream file_params;
	string filename_sequence;
	string filename_values;
	string filename_params;

	cout << "Введите имя файла, куда следует записать выборку: ";
	cin >> filename_sequence;
	cout << "Введите имя файла, куда следует записать плотности распределения: ";
	cin >> filename_values;
	cout << "Введите имя файла, куда следует записать параметры распределения: ";
	cin >> filename_params;

	file_sequence.open(filename_sequence);
	file_values.open(filename_values);
	file_params.open("empirical_params.txt");

	for (const pair<double, double>& pair : pairs)
	{
		file_sequence << pair.first << endl;
		file_values << pair.second << endl;
	}

	file_params << n << k;

	file_sequence.close();
	file_values.close();
	file_params.close();

	cout << "Выборка сохранена в файл " + filename_sequence << endl;
	cout << "Значения плотности сохранены в " + filename_values << endl;
	cout << "Параметры распределения сохранены в файл " + filename_params << endl << endl;
}

void Empirical::load_from_file(ifstream& file)
{
	string filename;

	cout << "Введите имя файла с выборкой элементов эмпирического распределения: ";
	cin >> filename;

	file.open(filename);
	if (!file)
		throw runtime_error("Ошибка: не удалось открыть файл");

	double x;
	while (!file.eof())
	{
		file >> x;
		x_s.push_back(x);
	}
	file.close();

	n = x_s.size();
	k = (int)floor(log2(n)) + 1;
	f_s = generate_values();
}


double Empirical::f(const double x) const
{
	int k = (int)floor(log2((double)n)) + 1;
	double min_x = *min_element(begin(x_s), end(x_s));
	double max_x = *max_element(begin(x_s), end(x_s));
	double delta = (max_x - min_x) / (double)k;

	for (int i = 0; i < k; i++)
		if (min_x + delta * i <= x && x < min_x + delta * (i + 1))
		{
			int n_i = count_if(x_s.begin(), x_s.end(), [i, k, min_x, max_x, delta](double x) { return i == k - 1 ? min_x + delta * (double)i <= x && x <= min_x + delta * (double)(i + 1) : min_x + delta * (double)i <= x && x < min_x + delta * (double)(i + 1); });
			return n_i / (n * delta);
		}

	return 0.0;
}

double Empirical::expected_value() const
{
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += x_s[i];
	}
	return sum / (double)n;
}

double Empirical::variance() const
{
	const double expected_val = expected_value();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_s[i] - expected_val, 2);

	return sum / n;
}

double Empirical::asymmetry() const
{
	const double expected_val = expected_value();
	const double variance_val = variance();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_s[i] - expected_val, 3);

	return sum / (n * pow(variance_val, 3 / 2));
}

double Empirical::kurtosis() const
{
	const double expected_val = expected_value();
	const double variance_val = variance();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		sum += pow(x_s[i] - expected_val, 4);
	}
	return (sum / (n * pow(variance_val, 2))) - 3;
}