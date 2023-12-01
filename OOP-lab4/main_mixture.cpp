#include "main_mixture.h"

void final_mixture(Mixture<Huber, Huber>* MX)
{
	ofstream x_selection;
	ofstream y_selection;
	ofstream params;
	vector<pair<double, double>> table;
	char option;
	int n;
	double x = 0;

	cout << "�������� �����:" << endl;
	cout << "1. ������� ��������� ����� ������������� �� �����" << endl;
	cout << "2. ��������� �������� ��������� ����� ������������� � ������������ �����" << endl;
	cout << "3. �������� ������� ��� �������" << endl;
	cout << "4. ������� ��������� ����� ������������� � ����" << endl;
	cout << "5. �����" << endl << endl;

	cin >> option;
	cout << endl;

	switch (option)
	{
	case '1':
		cout << "P = " << MX->get_p() << endl;
		cout << "V1 = " << MX->get_component1()->get_v() << endl;
		cout << "K1 = " << MX->get_component1()->get_k() << endl;
		cout << "Scale1 = " << MX->get_component1()->get_scale() << endl;
		cout << "Shift1 = " << MX->get_component1()->get_shift() << endl;
		cout << "V2 = " << MX->get_component2()->get_v() << endl;
		cout << "K2 = " << MX->get_component2()->get_k() << endl;
		cout << "Scale2 = " << MX->get_component2()->get_scale() << endl;
		cout << "Shift2 = " << MX->get_component2()->get_shift() << endl;
		cout << "M[X] = " << MX->M_Ksi() << endl;
		cout << "D[X] = " << MX->D_Ksi() << endl;
		cout << "As[X] = " << MX->asymmetry() << endl;
		cout << "Kurt[X] = " << MX->kurtosis() << endl;
		cout << "density(0) = " << MX->density(x) << endl << endl << endl;

		break;

	case '2':
		cout << "������� x: ";
		cin >> x;
		cout << endl << "density(" << x << ") = " << MX->density(x) << endl << endl;
		break;

	case '3':
		x_selection.open("selection1.txt");
		y_selection.open("selection2.txt");
		cout << "������� ����������� ������� n: ";
		cin >> n;
		table = MX->generate_pair(n);
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
		MX->save_to_file(params);
		cout << endl << endl;
		break;

	case '5':
		cout << endl << endl;
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl << endl;
		break;
	}

}

void next_mixture()
{
	char option;
	double p, v1, scale1, shift1, v2, scale2, shift2;
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
		MX = new Mixture<Huber, Huber>(new Huber(), new Huber(), 0.5);
		final_mixture(MX);
		break;

	case '2':
		cout << "������� p: " << endl; 
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
		final_mixture(MX);
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
		final_mixture(MX);
		break;

	case '4':
		break;

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}