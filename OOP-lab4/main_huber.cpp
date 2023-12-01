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

	cout << "��� ������ ������?" << endl;
	cout << "1. ������� ��������� ������������� �� �����" << endl;
	cout << "2. ��������� �������� ��������� ������������� � ������������ �����" << endl;
	cout << "3. �������� ������� ��� �������" << endl;
	cout << "4. ������� ��������� ������������� � ����" << endl;
	cout << "5. �����" << endl << endl;


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
		cout << "������� x: ";
		cin >> x;
		cout << endl << "density(" << x << ") = " << HB->density(x) << endl << endl;
		break;

	case '3':
		x_selection.open("selection1.txt");
		y_selection.open("selection2.txt");
		cout << "������� ����������� ������� n: ";
		cin >> n;
		table = HB->generate_pair(n);
		for (const pair<double, double>& pr : table)
		{
			x_selection << pr.first << endl;
			y_selection << pr.second << endl;
		}

		x_selection.close();
		x_selection.close();

		cout << endl << "�������� ������� �������� � ����� selection1.txt � selection2.txt" << endl << endl;
		break;

	case '4':
		HB->save_to_file(params);
		cout << endl << endl;
		break;

	case '5':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}

void next_huber()
{
	char option;
	double v, scale, shift;
	Huber* HB;
	ifstream file;

	cout << "�������� �����: " << endl;
	cout << "1. ����������� ������������� �������" << endl;
	cout << "2. ������ ������������ ��������� � ����������" << endl;
	cout << "3. ������� ��������� �� �����" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "������� �������� v: ";
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
		cout << "������� v: " << endl;
		cin >> v;
		cout << "������� scale: " << endl;
		cin >> scale;
		cout << "������� shift: " << endl;
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
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}













