#include "function.h"
#include "Distribution.h"

//--------------------------------------------------------------------------huber_distribution--------------------------------------------------------------------------------//
HuberD::HuberD(double v, double scale, double shift)
{
	if (v <= 0 || scale <= 0)
		throw invalid_argument("Ошибка: один или несколько параметров не корректны");

	this->v = v;
	this->k = K(v);
	this->scale = scale;
	this->shift = shift;
}

HuberD::HuberD(ifstream& file)
{
	load_file(file);
}

double HuberD::get_v() const
{
	return this->v;
}

double HuberD::get_k() const
{
	return this->k;
}

double HuberD::get_scale() const
{
	return this->scale;
}

double HuberD::get_shift() const
{
	return this->shift;
}

void HuberD::set_v(const double v)
{
	if (v <= 0)
		throw invalid_argument("Ошибка: v <= 0");

	this->v = v;
	this->k = K(this->v);
}

void HuberD::set_scale(const double scale)
{
	if (scale <= 0)
		throw invalid_argument("Ошибка: λ <= 0");
	this->scale = scale;
}

void HuberD::set_shift(const double shift)
{
	this->shift = shift;
}

double HuberD::Huber(const double x) const
{
	if (abs((x - this->shift) / this->scale) <= this->v)
	{
		return (1. / (sqrt(2. * M_PI) * this->k) * exp(-pow((x - this->shift) / this->scale, 2.) / 2.)) / this->scale;
	}
	if (abs((x - this->shift) / this->scale) > this->v)
	{
		return (1. / (sqrt(2. * M_PI) * this->k) * exp(pow(this->v, 2.) / 2. - this->v * abs((x - this->shift) / this->scale))) / this->scale;
	}
}
double HuberD::phi(double x) const//Ф(x)
{
	return 0.5 * (1. + erf(x / sqrt(2.)));
}
double HuberD::phi_lower(double x) const//ф(x)
{
	return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}
double HuberD::Mksi_huber() const//мат ожидание
{
	return this->shift;
}
double HuberD::Dksi_huber() const//дисперсия
{
	return 1. + 2. * phi_lower(this->v) * (pow(this->v, 2.) + 2.) / (pow(this->v, 3.) * this->k);
}
double HuberD::asymmetry_huber() const//ассиметрия
{
	return 0.;
}
double HuberD::kurtosis_huber() const//коэфф эксцесса
{
	return (3. * (2. * phi(this->v) - 1.) + 2. * phi_lower(this->v) * (24. / pow(this->v, 5.) + 24. / pow(this->v, 3.) + 12. /
		this->v + this->v)) / (pow(Dksi_huber(), 2.) * this->k) - 3.;
}
double HuberD::P() const//вероятности попадания в центральный интервал 
{
	return (2. * phi(this->v) - 1.) / this->k;
}
double HuberD::K(double v) const//значение K зная V
{
	return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}
//-----------------------------------------------------//
double HuberD::algorithm() const
{
	std::random_device rd;
	std::default_random_engine gen(rd());
	std::uniform_real_distribution<> d(0, 1);
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

void HuberD::load_file(ifstream& file)
{
	string filename;
	double v, scale, shift;

	file.open("input.txt");

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

void HuberD::save_file(ofstream& file) const
{
	string filename = "output.txt";

	file.open(filename);
	file << this->v << endl << this->scale << endl << this->shift;
	file.close();
}

vector<double> HuberD::selection(const int n) const
{
	vector<double> res;

	for (int i = 0; i < n; i++)
	{
		double x = algorithm();
		res.push_back(x);
	}

	sort(res.begin(), res.end());

	return res;
}

vector<pair<double, double>> HuberD::generate_pair(const int n, const vector<double>& x_selection) const
{
	vector<pair<double, double>> res;
	vector<double> sequence;

	if (x_selection.empty())//если строка пуста заполни ее
		sequence = selection(n);
	else
		sequence = x_selection;

	for (const double& x : sequence)
	{
		double y = Huber(x);
		res.push_back(make_pair(x, y));
	}

	return res;
}





//--------------------------------------------------------------------------mixture_distribution--------------------------------------------------------------------------------//

template<class Distribution1, class Distribution2>
class Mixture : public IDistribution, public IPersistent
{
private:

	Distribution1* HD1;

	Distribution2* HD2;

	double p;

public:

	Mixture(Distribution1* D1, Distribution2* D2, double p) : HD1(HD1), HD2(HD2), p(p) {};

	void load_file(ifstream& file) override;

	Distribution1* get_component1();//+

	Distribution2* get_component2();//+

	double get_p() const;//+

	void set_p(const double p);//+

	void save_to_file(ofstream& file) override;

	double H_Mixture(const double x) const;//+

	double Mksi_mixture() const;//+

	double Dksi_Mixture() const;//+

	double asymmetry_mixture() const;//+

	double kurtosis_mixture() const;//+

	double algorithm_mixture() const;//+

	vector<double> selection_mixture(const int n) const;

	vector<pair<double, double>> generate_pair_mixture(const int n, const vector<double>& x_selection = {}) const;


};


//Mixture::Mixture(HuberD* _HB1, HuberD* _HB2, double _p) :
//	HD1(_HB1), HD2(_HB2), p(_p > 0 and _p < 1 ? _p : throw invalid_argument("Ошибка: один или несколько параметров некорректны")) {}


template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::load_file(ifstream& file)
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
Distribution1* Mixture<Distribution1, Distribution2>::get_component1()
{
	return HD1;
}

template <class Distribution1, class Distribution2>
Distribution2* Mixture<Distribution1, Distribution2>::get_component2()
{
	return HD2;
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::get_p() const
{
	return p;
}

template <class Distribution1, class Distribution2>
void Mixture<Distribution1, Distribution2>::set_p(const double p)
{
	if (p > 1 || p < 0)
		throw invalid_argument("Ошибка: один параметро некорректен");

	this->p = p;
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::H_Mixture(const double x) const//плотность
{
	return (1 - p) * HD1->Huber(x) + p * HD2->Huber(x);
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::Mksi_mixture() const
{
	return (1 - p) * HD1->Mksi_huber() + p * HD2->Mksi_huber();
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::Dksi_Mixture() const
{
	return (1 - p) * (pow(HD1->Mksi_huber(), 2) + HD1->Dksi_huber()) +
		p * (pow(HD2->Mksi_huber(), 2) + HD2->Dksi_huber()) -
		pow(Mksi_mixture(), 2);


}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::asymmetry_mixture() const
{
	return (1 / pow(Dksi_Mixture(), 3 / 2)) *
		((1 - p) *
			(pow(HD1->Mksi_huber() - Mksi_mixture(), 3) + 3 *
				(HD1->Mksi_huber() - Mksi_mixture()) * HD1->Dksi_huber() +
				pow(HD1->Dksi_huber(), 3 / 2) * HD1->asymmetry_huber()) +
			(p) *
			(pow(HD2->Mksi_huber() - Mksi_mixture(), 3) + 3 *
				(HD2->Mksi_huber() - Mksi_mixture()) * HD2->Dksi_huber() +
				pow(HD2->Dksi_huber(), 3 / 2) * HD2->asymmetry_huber()
				));
}

template <class Distribution1, class Distribution2>
double Mixture<Distribution1, Distribution2>::kurtosis_mixture() const
{
	return (1 / pow(Dksi_Mixture(), 2)) * (
		(1 - p) * (pow(HD1->Mksi_huber() - Mksi_mixture(), 4) +
			6 * pow(HD1->Mksi_huber() - Mksi_mixture(), 2) * HD1->Dksi_huber() +
			4 * (HD1->Mksi_huber() - Mksi_mixture()) * pow(HD1->Dksi_huber(), 3 / 2) * HD1->asymmetry_huber() +
			pow(HD1->Dksi_huber(), 2) * HD1->kurtosis_huber()) +
		(p) * (pow(HD2->Mksi_huber() - Mksi_mixture(), 4) +
			6 * pow(HD2->Mksi_huber() - Mksi_mixture(), 2) * HD2->Dksi_huber() +
			4 * (HD2->Mksi_huber() - Mksi_mixture()) * pow(HD2->Dksi_huber(), 3 / 2) * HD2->asymmetry_huber() +
			pow(HD2->Dksi_huber(), 2) * HD2->kurtosis_huber()) - 3);
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
double Mixture<Distribution1, Distribution2>::algorithm_mixture() const {//моделирует случайную величину для смеси
	random_device rd;
	default_random_engine gen(rd());
	uniform_real_distribution<> d(0, 1);

	double r = d(gen);

	if (r > p)
		return HD1->algorithm();
	else
		return HD2->algorithm();
}

template <class Distribution1, class Distribution2>
vector<double> Mixture<Distribution1, Distribution2>::selection_mixture(const int n) const//формирование выборки случайно величины по закону смеси
{
	vector<double> sequence;
	for (int i = 0; i < n; i++)
	{
		double x = algorithm_mixture();
		sequence.push_back(x);
	}

	sort(sequence.begin(), sequence.end());
	return sequence;
}


template <class Distribution1, class Distribution2>
vector<pair<double, double>> Mixture<Distribution1, Distribution2>::generate_pair_mixture(const int n, const vector<double>& x_selection) const
{
	vector<double> sequence;
	vector<pair<double, double>> table;

	if (x_selection.empty())
		sequence = selection_mixture(n);
	else
		sequence = x_selection;

	for (const double& x : sequence)
	{
		double y = H_Mixture(x);
		table.push_back(make_pair(x, y));
	}

	return table;
}

//--------------------------------------------------------------------------empirical_distribution--------------------------------------------------------------------------------//

Empirical::Empirical(const IDistribution* D, int _n, int _k) :
	n(_n > 1 ? _n : throw invalid_argument("Некорректный аргумент")), k(_k > 2 ? _k : (int)floor(log2(_n)) + 1)
{
	x_selection = D->generate_sequence(_n);
	f_selection = generate_values();
}

Empirical::Empirical(const Empirical* EM) : n(EM->n > 1 ? EM->n : throw invalid_argument("Некорректный аргумент")), k(EM->k > 2 ? EM->k : (int)floor(log2(EM->n) + 1)), x_selection(EM->x_selection), f_selection(EM->f_selection)
{
	f_selection = generate_values();
}

Empirical::Empirical(const int _n, const int _k) : n(_n > 1 ? _n : throw invalid_argument("Некорректный аргумент")), k(_k > 2 ? _k : (int)floor(log2(_n) + 1))
{
	x_selection = generate_sequence(n);
	f_selection = generate_values();
}

Empirical::Empirical(const vector<double>& x_s) : n(x_s.size()), k((int)floor(log2(x_s.size())) + 1)
{
	this->x_selection = x_s;
	f_selection = generate_values();
}

Empirical::Empirical(ifstream& file)
{
	load_file(file);
}

Empirical::~Empirical()
{
	x_selection.clear();
	f_selection.clear();
}

Empirical& Empirical::operator=(const Empirical& EM)
{
	if (this == &EM)
		return *this;

	x_selection = EM.x_selection;
	f_selection = EM.f_selection;
	n = EM.n;
	k = EM.k;
	return *this;
}

double Empirical::algorithm_empirical() const
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
		result.push_back(algorithm_empirical());

	sort(result.begin(), result.end());

	return result;
}


vector<double> Empirical::generate_values() const
{
	vector<double> result;

	for (const double& x : x_selection)
		result.push_back(H_Empirical(x));

	return result;
}

vector<pair<double, double>> Empirical::generate_pair(const int n, const vector<double>& x_selection) const
{
	vector<pair<double, double >> result;

	for (int i = 0; i < n; i++)
		result.push_back(make_pair(x_selection[i], f_selection[i]));

	return result;
}

vector<double> Empirical::get_x_selection() const
{
	return x_selection;
}

vector<double> Empirical::get_f_selection() const
{
	return f_selection;
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
	vector<pair<double, double>> pairs = generate_pair(n);

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

void Empirical::load_file(ifstream& file)
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
		x_selection.push_back(x);
	}
	file.close();

	n = x_selection.size();
	k = (int)floor(log2(n)) + 1;
	f_selection = generate_values();
}


double Empirical::H_Empirical(const double x) const
{
	int k = (int)floor(log2((double)n)) + 1;
	double min_x = *min_element(begin(x_selection), end(x_selection));
	double max_x = *max_element(begin(x_selection), end(x_selection));
	double delta = (max_x - min_x) / (double)k;

	for (int i = 0; i < k; i++)
		if (min_x + delta * i <= x && x < min_x + delta * (i + 1))
		{
			int n_i = count_if(x_selection.begin(), x_selection.end(), [i, k, min_x, max_x, delta](double x) { return i == k - 1 ? min_x + delta * (double)i <= x && x <= min_x + delta * (double)(i + 1) : min_x + delta * (double)i <= x && x < min_x + delta * (double)(i + 1); });
			return n_i / (n * delta);
		}

	return 0.0;
}

double Empirical::Mn() const
{
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += x_selection[i];
	}
	return sum / (double)n;
}


double Empirical::Dn() const//+
{
	const double Mn_ksi = Mn();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		double diff = x_selection[i] - Mn_ksi;
		sum += diff * diff;
	}

	return sum / n;
}

double Empirical::asymmetry_empirical() const
{
	const double Mn_ksi = Mn();
	const double Dn_ksi = Dn();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_selection[i] - Mn_ksi, 3);

	return sum / (n * pow(Dn_ksi, 3 / 2));
}

double Empirical::kurtosis_empirical() const
{
	const double Mn_ksi = Mn();
	const double Dn_ksi = Dn();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		sum += pow(x_selection[i] - Mn_ksi, 4);
	}
	return (sum / (n * pow(Dn_ksi, 2))) - 3;
}

