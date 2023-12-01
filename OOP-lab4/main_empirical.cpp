#include "main_empirical.h"
#include "mixture.cpp"
#include "function.h"
#include "lib.h"

void final_empirical(Empirical* EM)
{
	ofstream file;
	char option;
	double x = 0;

	cout << "�������� �����:" << endl;
	cout << "1. ������� ��������� ������������� ������������� �� �����" << endl;
	cout << "2. ��������� �������� ��������� ����� ������������� � ������������ �����" << endl;
	cout << "3. ������� ������� ��������� � ��������� ������������� ������������� � ����" << endl;
	cout << "4. �����" << endl << endl;

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
		cout << "������� x: ";
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
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "�������� �����: " << endl;
	cout << "1. ����������� ������������� �������" << endl;
	cout << "2. ������ ������������ ��������� � ����������" << endl;
	cout << "3. ������� ��������� �� �����" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(HB, n, k);
		final_empirical(EM);

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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;

		EM = new Empirical(HB, n, k);
		final_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "�������� �����: " << endl;
	cout << "1. ����������� ��������� ����� �������������" << endl;
	cout << "2. ���� � ����������" << endl;
	cout << "3. ���� �� �����" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		MX = new Mixture<Huber, Huber>(new Huber(), new Huber, 0.5);

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MX, n, k);

		final_empirical(EM);
		break;

	case '2':
		cout << "�������: " << endl;
		cin >> p;
		cout << "������� v1: " << endl;
		cin >> v1;
		cout << "������� scale1: " << endl;
		cin >> scale1;
		cout << "������� shift1: " << endl;
		cin >> shift1;
		cout << "������� v2: " << endl;
		cin >> v2;
		cout << "������� scale2: " << endl;
		cin >> scale2;
		cout << "������� shift2: " << endl;
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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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

		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
		cin >> k;
		cout << endl << endl;

		EM = new Empirical(MX, n, k);
		final_empirical(EM);
		break;

	case '4':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "1. ������ ��������� ������������� �� �����" << endl;
	cout << "2. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		cout << "������� ��� ����� � ����������� �������������: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("������: �� ������� ������� ����");

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
		cerr << endl << "������: ������������ ��������" << endl << endl;
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

	cout << "1. ������ ������� ��������� �� �����" << endl;
	cout << "2. �����" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "������� ��� ����� � ����������� �������������: ";
		cin >> filename;
		cout << endl << endl;

		file.open(filename);

		if (!file)
			throw runtime_error("������: �� ������� ������� ����");

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
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;

	}

}

void empirical()
{
	Empirical* EM;
	ifstream file;
	int n, k;
	char option;

	cout << "������ ����� ��������� ��� ������������� �������������?" << endl;
	cout << "1. �������� �������������" << endl;
	cout << "2. ����� �������������" << endl;
	cout << "3. �� ���� ������������� ������������� �������������" << endl;
	cout << "4. �� ���� ������������ ������� ��������� ��������� ��������" << endl;
	cout << "5. ����������� ���������" << endl;
	cout << "6. ������ ��������� �� �����" << endl;
	cout << "7. �����" << endl << endl;

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
		cout << "������� n: ";
		cin >> n;
		cout << "������� k (0, ���� ������ ���������� k �� ������� ����������): ";
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
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}