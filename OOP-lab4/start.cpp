#include "main_huber.h"
#include "main_empirical.h"
#include "main_mixture.h"
#include "lib.h"


void step1()
{
	char option;

	cout << "�������������: " << endl;
	cout << "1. �������� �������������" << endl;
	cout << "2. ����� �������������" << endl;
	cout << "3. ������������ �������������" << endl;
	cout << "4. �����" << endl << endl;

	cin >> option;
	cout << endl << endl;

	switch (option)
	{
	case '1':
		next_huber();
		break;

	case '2':
		next_mixture();
		break;

	case '3':
		empirical();
		break;

	case '4':
		cout << "���������� ������" << endl << endl;
		exit(0);

	default:
		cerr << endl << "������: ������������ ��������" << endl << endl;
		break;
	}
}

void start()
{
	while (true)
		step1();
}