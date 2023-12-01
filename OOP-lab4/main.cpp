#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
//#include "Interface.h"

/*void demo()
{
	cout << "7. Тест позднего связывания:" << endl << endl;

	Primary* HB1 = new Primary(1, 1, 2);
	Primary* HB2 = new Primary(1, 2, 6);
	Primary* HB3 = new Primary(1, 2, -6);
	Primary* HB4 = new Primary(1, 1, -2);

	Mixture<Primary, Primary>* MX1 = new Mixture<Primary, Primary>(HB1, HB2, 0.5);
	Mixture<Primary, Primary>* MX2 = new Mixture<Primary, Primary>(HB3, HB4, 0.5);
	Mixture<Mixture<Primary, Primary>, Mixture<Primary, Primary>>* MX = new Mixture<Mixture<Primary, Primary>, Mixture<Primary, Primary>>(MX1, MX2, 0.5);


	cout << "Параметр формы 1: " << MX->get_component1()->get_component1()->get_v() << endl;
	cout << "Параметр масштаба 1: " << MX->get_component1()->get_component1()->get_scale() << endl;
	cout << "Параметр сдвига 1: " << MX->get_component1()->get_component1()->get_shift() << endl << endl;

	cout << "Параметр формы 2: " << MX->get_component1()->get_component2()->get_v() << endl;
	cout << "Параметр масштаба 2: " << MX->get_component1()->get_component2()->get_scale() << endl;
	cout << "Параметр сдвига 2: " << MX->get_component1()->get_component2()->get_shift() << endl << endl;

	cout << "Параметр формы 3: " << MX->get_component2()->get_component1()->get_v() << endl;
	cout << "Параметр масштаба 3: " << MX->get_component2()->get_component1()->get_scale() << endl;
	cout << "Параметр сдвига 3: " << MX->get_component2()->get_component1()->get_shift() << endl << endl;

	cout << "Параметр формы 4: " << MX->get_component2()->get_component2()->get_v() << endl;
	cout << "Параметр масштаба 4: " << MX->get_component2()->get_component2()->get_scale() << endl;
	cout << "Параметр сдвига 4: " << MX->get_component2()->get_component2()->get_shift() << endl << endl;

	cout << "Теоретические характеристики:" << endl;
	cout << "Математическое ожидание: " << MX->expected_value() << endl;
	cout << "Дисперсия: " << MX->variance() << endl;
	cout << "Коэффиицент асимметрии: " << MX->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX->kurtosis() << endl << endl;

	Empirical* EM1 = new Empirical(MX, 10000, 0);
	Empirical* EM2 = new Empirical(MX, 10000, 0);

	cout << "Эмпирические характеристики:" << endl;
	cout << "Математическое ожидание: " << MX->expected_value() << endl;
	cout << "Дисперсия: " << MX->variance() << endl;
	cout << "Коэффициент асимметрии: " << MX->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX->kurtosis() << endl << endl;


	vector<pair<double, double>> s1 = MX->generate_table_of_values(10000, MX->generate_sequence(10000));
	vector<pair<double, double>> s2 = EM1->generate_table_of_values(10000, EM1->generate_sequence(10000));
	vector<pair<double, double>> s3 = EM2->generate_table_of_values(10000, EM2->generate_sequence(10000));

	ofstream file1;
	ofstream file2;
	ofstream file3;

	file1.open("theoretical.txt");
	file2.open("empirical.txt");
	file3.open("empirical2.txt");

	for (auto& s : s1)
		file1 << s.second << endl;
	for (auto& s : s2)
		file2 << s.second << endl;
	for (auto& s : s3)
		file3 << s.second << endl;

	file1.close();
	file2.close();
	file3.close();
}*/

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");

	//demo();

	//start();

	int result = Catch::Session().run(argc, argv);
	return result;

	return 0;
}