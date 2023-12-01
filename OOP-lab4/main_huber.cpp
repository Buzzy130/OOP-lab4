#include "main_huber.h"

void final_primary(Huber* HB)
{
	ofstream x_selection;
	ofstream y_selection;
	ofstream params;
	vector<pair<double, double>> table;
	char option;
	int n;
	double x = 0;

	cout << "Что делаем дальше?" << endl;
	cout << "1. Вывести параметры распределения на экран" << endl;
	cout << "2. Вычислить значение плотности распределения в произвольной точке" << endl;
	cout << "3. Получить выборку для анализа" << endl;
	cout << "4. Вывести параметры распределения в файл" << endl;
	cout << "5. Назад" << endl << endl;


	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "v = " << HB->get_v() << endl;
		cout << "K = " << HB->get_k() << endl;
		cout << "Scale = " << HB->get_scale() << endl;
		cout << "Shift = " << HB->get_shift() << endl;
		cout << "M[X] = " << HB->M_Ksi() << endl;
		cout << "D[X] = " << HB->D_Ksi() << endl;
		cout << "As[X] = " << HB->asymmetry() << endl;
		cout << "Kurt[X] = " << HB->kurtosis() << endl;
		cout << "P = " << HB->P() << endl;
		cout << "density(0) = " << HB->density(0) << endl << endl << endl;

		break;

	case '2':
		cout << "Введите x: ";
		cin >> x;
		cout << endl << "density(" << x << ") = " << HB->density(x) << endl << endl;
		break;

	case '3':
		x_selection.open("selection1.txt");
		y_selection.open("selection2.txt");
		cout << "Введите размерность выборки n: ";
		cin >> n;
		table = HB->generate_pair(n);
		for (const pair<double, double>& pr : table)
		{
			x_selection << pr.first << endl;
			y_selection << pr.second << endl;
		}

		x_selection.close();
		x_selection.close();

		cout << endl << "Значения выборки записаны в файлы selection1.txt и selection2.txt" << endl << endl;
		break;

	case '4':
		HB->save_to_file(params);
		cout << endl << endl;
		break;

	case '5':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void next_huber()
{
	char option;
	double v, scale, shift;
	Huber* HB;
	ifstream file;

	cout << "Выберите опцию: " << endl;
	cout << "1. Стандартное распределение Хьюбера" << endl;
	cout << "2. Ввести произвольные параметры с клавиатуры" << endl;
	cout << "3. Считать параметры из файла" << endl;
	cout << "4. Назад" << endl << endl;

	cin >> option;
	cout << endl;

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

		final_primary(HB);
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

		final_primary(HB);
		break;

	case '3':
		try
		{
			HB = new Huber(file);
		}
		catch (exception e)
		{
			cerr << endl << e.what() << endl << endl;
			break;
		}

		final_primary(HB);
		break;

	case '4':
		break;

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}













