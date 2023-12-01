#include "main_empirical.h"
#include "mixture.cpp"
#include "function.h"
#include "lib.h"

void final_empirical(Empirical* EM)
{
	ofstream file;
	char option;
	double x = 0;

	cout << "Выберите опцию:" << endl;
	cout << "1. Вывести параметры эмпирического распределения на экран" << endl;
	cout << "2. Вычислить значение плотности смеси распределений в произвольной точке" << endl;
	cout << "3. Вывести выборку элементов и параметры эмпирического распределения в файл" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "n = " << EM->get_n() << endl;
		cout << "K = " << EM->get_k() << endl;
		cout << "M[X] = " << EM->M_Ksi() << endl;
		cout << "D[X] = " << EM->D_Ksi() << endl;
		cout << "As[X] = " << EM->asymmetry() << endl;
		cout << "Kurt[X] = " << EM->kurtosis() << endl;
		cout << "density(0) = " << EM->density(x) << endl << endl << endl;

		break;

	case '2':
		cout << "Введите x: ";
		cin >> x;
		cout << endl << "density(" << x << ") = " << EM->density(x) << endl << endl;
		break;

	case '3':
		EM->save_to_file(file);
		cout << endl << endl;
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}


void empirical_primary()
{
	Empirical* EM;
	char option;
	double v, scale, shift;
	int n, k;
	Huber* HB;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартное распределение Хьюбера" << endl;
	cout << "2. Ввести произвольные параметры с клавиатуры" << endl;
	cout << "3. Считать параметры из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "Введите параметр v: ";
		cin >> v;
		cout << endl << endl;
		try
		{
			HB = new Huber(v);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(HB, n, k);
		final_empirical(EM);

		break;

	case '2':
		cout << "Введите v: " << endl;
		cin >> v;
		cout << "Введите scale: " << endl;
		cin >> scale;
		cout << "Введите shift: " << endl;
		cin >> shift;
		cout << endl << endl;
		try
		{
			HB = new Huber(v, scale, shift);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;
		EM = new Empirical(HB, n, k);
		final_empirical(EM);

		break;

	case '3':
		try
		{
			HB = new Huber(file);
		}
		catch (exception e)
		{
			cerr << endl << e.what() << endl;
			break;
		}

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;

		EM = new Empirical(HB, n, k);
		final_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void empirical_mixture()
{
	Empirical* EM;
	char option;
	double p, v1, scale1, shift1, v2, scale2, shift2;
	int n, k;
	Mixture<Huber, Huber>* MX;
	Huber* HB1;
	Huber* HB2;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартные параметры смеси распределений" << endl;
	cout << "2. Ввод с клавиатуры" << endl;
	cout << "3. Ввод из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		MX = new Mixture<Huber, Huber>(new Huber(), new Huber, 0.5);

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MX, n, k);

		final_empirical(EM);
		break;

	case '2':
		cout << "Введите: " << endl;
		cin >> p;
		cout << "Введите v1: " << endl;
		cin >> v1;
		cout << "Введите scale1: " << endl;
		cin >> scale1;
		cout << "Введите shift1: " << endl;
		cin >> shift1;
		cout << "Введите v2: " << endl;
		cin >> v2;
		cout << "Введите scale2: " << endl;
		cin >> scale2;
		cout << "Введите shift2: " << endl;
		cin >> shift2;
		try
		{
			HB1 = new Huber(v1, scale1, shift1);
			HB2 = new Huber(v2, scale2, shift2);
			MX = new Mixture<Huber, Huber>(HB1, HB2, p);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MX, n, k);
		final_empirical(EM);
		break;

	case '3':
		MX = new Mixture<Huber, Huber>(new Huber(), new Huber(), 0.5);
		try
		{
			MX->load_from_file(file);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MX, n, k);
		final_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void empirical_empirical()
{
	Empirical* EM;
	char option;
	ifstream file;
	string filename;
	int n, k;

	cout << "1. Ввести параметры распределения из файла" << endl;
	cout << "2. Назад" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "Введите имя файла с параметрами распределения: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("Ошибка: не удалось открыть файл");

		while (!file.eof())
			file >> n >> k;

		file.close();

		try
		{
			EM = new Empirical(n, k);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		cout << endl << endl;

		final_empirical(EM);
		break;

	case '2':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;

	}
}

void empirical_sequence()
{
	Empirical* EM;
	char option;
	ifstream file;
	string filename;
	double x;
	vector<double> x_s;

	cout << "1. Ввести выборку элементов из файла" << endl;
	cout << "2. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "Введите имя файла с параметрами распределения: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("Ошибка: не удалось открыть файл");

		while (!file.eof())
		{
			file >> x;
			x_s.push_back(x);
		}

		file.close();

		try
		{
			EM = new Empirical(x_s);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}

		cout << endl << endl;
		final_empirical(EM);
		break;

	case '2':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;

	}

}

void empirical()
{
	Empirical* EM;
	ifstream file;
	int n, k;
	char option;

	cout << "Откуда берем параметры для эмпирического распределения?" << endl;
	cout << "1. Основное распределение" << endl;
	cout << "2. Смесь распределений" << endl;
	cout << "3. На базе существующего эмпирического распределения" << endl;
	cout << "4. На базе существующей выборки элементов случайной величины" << endl;
	cout << "5. Собственные параметры" << endl;
	cout << "6. Ввести параметры из файла" << endl;
	cout << "7. Назад" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		empirical_primary();
		break;

	case '2':
		empirical_mixture();
		break;

	case '3':
		empirical_empirical();
		break;

	case '4':
		empirical_sequence();
		break;

	case '5':
		cout << "Введите n: ";
		cin >> n;
		cout << "Введите k (0, если хотите рассчитать k по формуле Стерджесса): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(n, k);
		final_empirical(EM);
		break;

	case '6':
		try {
			EM = new Empirical(file);
		}
		catch (exception e)
		{
			cerr << e.what() << endl << endl;
			break;
		}
		final_empirical(EM);
		break;

	case '7':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}