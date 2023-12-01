#include "main_huber.h"
#include "main_empirical.h"
#include "main_mixture.h"
#include "lib.h"


void step1()
{
	char option;

	cout << "Распределение: " << endl;
	cout << "1. Основное распределение" << endl;
	cout << "2. Смесь распределений" << endl;
	cout << "3. Эмпирическое распределение" << endl;
	cout << "4. Выход" << endl << endl;

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
		cout << "Завершение работы" << endl << endl;
		exit(0);

	default:
		cerr << endl << "Ошибка: Некорректный параметр" << endl << endl;
		break;
	}
}

void start()
{
	while (true)
		step1();
}