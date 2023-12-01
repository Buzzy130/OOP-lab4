#include "function.h"

Huber::Huber(double v, double scale, double shift)
{
    if (v <= 0 || scale <= 0)
        throw invalid_argument("Ошибка: один или несколько параметров не корректны");

    this->v = v;
    this->k = K(v);
    this->scale = scale;
    this->shift = shift;
}

Huber::Huber(ifstream& file)
{
    load_from_file(file);
}

void Huber::save_to_file(ofstream& file)
{
    string filename;
    cout << "Введите имя файла, куда следует записать параметры основного распределения: ";
    cin >> filename;

    file.open(filename);
    file << this->v << endl << this->scale << endl << this->shift;
    file.close();

    cout << "Параметры основного распределения сохранены в файл " + filename << endl;
}

void Huber::load_from_file(ifstream& file)
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


double Huber::get_v() const
{
    return this->v;
}

double Huber::get_k() const
{
    return this->k;
}


void Huber::set_v(double v)
{
    if (v <= 0)
        throw invalid_argument("Ошибка: v <= 0");

    this->v = v;
    this->k = K(this->v);
}

double Huber::get_scale() const
{
    return this->scale;
}

void Huber::set_scale(const double scale)
{
    if (scale <= 0)
        throw invalid_argument("Ошибка: scale <= 0");
    this->scale = scale;
}

double Huber::get_shift() const
{
    return this->shift;
}

void Huber::set_shift(const double shift)
{
    this->shift = shift;
}

double Huber::density(const double x) const
{
	double z = (x - this->shift) / this->scale;
	double abs_z = abs(z);

	double numerator = (abs_z <= this->v) ? exp(-pow(z, 2.) / 2.) : exp(pow(this->v, 2.) / 2. - this->v * abs_z);
	double denominator = sqrt(2. * M_PI) * this->k * this->scale;

	return numerator / denominator;
}

double Huber::phi(double x) const
{
    return 0.5 * (1. + erf(x / sqrt(2.)));
}

double Huber::phi_lower(double x) const
{
    return 1. / sqrt(2. * M_PI) * exp(-1. / 2. * pow(x, 2.));
}

double Huber::M_Ksi() const
{
    return this->shift;
}

double Huber::D_Ksi() const
{
    return 1. + 2. * phi_lower(this->v) * (pow(this->v, 2.) + 2.) / (pow(this->v, 3.) * this->k);
}

double Huber::asymmetry() const
{
    return 0.;
}

double Huber::kurtosis() const
{
    return (3. * (2. * phi(this->v) - 1.) + 2. * phi_lower(this->v) * (24. / pow(this->v, 5.) + 24. / pow(this->v, 3.) + 12. / this->v + this->v)) / (pow(D_Ksi(), 2.) * this->k) - 3.;
}

double Huber::P() const
{
    return (2. * phi(this->v) - 1.) / this->k;
}

double Huber::K(const double v) const
{
    return 2. / v * phi_lower(v) + 2. * phi(v) - 1.;
}

double Huber::algorithm() const
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

vector<double> Huber::selection(const int n) const
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

vector<pair<double, double>> Huber::generate_pair(const int n, const vector<double>& x_s) const
{
    vector<pair<double, double>> res;
    vector<double> sequence;

    if (x_s.empty())
        sequence = selection(n);
    else
        sequence = x_s;

    for (const double& x : sequence)
    {
        double y = density(x);
        res.push_back(make_pair(x, y));
    }

    return res;
}
//===============================================================================//

Empirical::Empirical(const IDistribution* D, int _n, int _k)
{
	if (_n <= 1)
		throw invalid_argument("Некорректный аргумент");

	n = _n;
	k = (_k > 2) ? _k : (int)floor(log2(_n)) + 1;

	x_selection = D->selection(n);
	f_selection = generate_f_selection();
}

Empirical::Empirical(const Empirical* EM)
{
	if (EM->n <= 1)
		throw invalid_argument("Некорректный аргумент");

	n = EM->n;
	k = (EM->k > 2) ? EM->k : (int)floor(log2(EM->n) + 1);
	x_selection = EM->x_selection;
	f_selection = generate_f_selection();
}

Empirical::Empirical(const int _n, const int _k)
{
	if (_n <= 1)
		throw invalid_argument("Некорректный аргумент");

	n = _n;
	k = (_k > 2) ? _k : (int)floor(log2(_n) + 1);

	x_selection = selection(n);
	f_selection = generate_f_selection();
}

Empirical::Empirical(const vector<double>& x_selection)
{
	n = x_selection.size();
	k = (int)floor(log2(n) + 1);
	this->x_selection = x_selection;
	f_selection = generate_f_selection();
}

Empirical::Empirical(ifstream& file)
{
	load_from_file(file);
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

double Empirical::algorithm() const
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

vector<double> Empirical::selection(const int n) const
{
	vector<double> result;

	for (int i = 0; i < n; i++)
		result.push_back(algorithm());

	sort(result.begin(), result.end());

	return result;
}


vector<double> Empirical::generate_f_selection() const
{
	vector<double> result;

	for (const double& x : x_selection)
		result.push_back(density(x));

	return result;
}

vector<pair<double, double>> Empirical::generate_pair(const int n, const vector<double>& x_s) const
{
	vector<pair<double, double >> result;

	for (int i = 0; i < n; i++)
		result.push_back(make_pair(x_s[i], f_selection[i]));

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
		x_selection.push_back(x);
	}
	file.close();

	n = x_selection.size();
	k = (int)floor(log2(n)) + 1;
	f_selection = generate_f_selection();
}


double Empirical::density(const double x) const
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

double Empirical::M_Ksi() const
{
	double sum = 0;
	for (int i = 0; i < n; ++i) {
		sum += x_selection[i];
	}
	return sum / (double)n;
}

double Empirical::D_Ksi() const
{
	const double expected_val = M_Ksi();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_selection[i] - expected_val, 2);

	return sum / n;
}

double Empirical::asymmetry() const
{
	const double expected_val = M_Ksi();
	const double variance_val = D_Ksi();
	double sum = 0;

	for (int i = 0; i < n; ++i)
		sum += pow(x_selection[i] - expected_val, 3);

	return sum / (n * pow(variance_val, 3 / 2));
}

double Empirical::kurtosis() const
{
	const double expected_val = M_Ksi();
	const double variance_val = D_Ksi();
	double sum = 0;

	for (int i = 0; i < n; ++i) {
		sum += pow(x_selection[i] - expected_val, 4);
	}
	return (sum / (n * pow(variance_val, 2))) - 3;
}