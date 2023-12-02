#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "start.h"

/*void demo()
{
	cout << "7. Тест позднего связывания:" << endl << endl;

	Huber* HB1 = new Huber(1, 1, 0);
	Huber* HB2 = new Huber(1, 1, 0);
	Huber* HB3 = new Huber(1, 1, 0);
	Huber* HB4 = new Huber(1, 1, 0);

	Mixture<Huber, Huber>* MX1 = new Mixture<Huber, Huber>(HB1, HB2, 0.5);
	Mixture<Huber, Huber>* MX2 = new Mixture<Huber, Huber>(HB3, HB4, 0.5);
	Mixture<Mixture<Huber, Huber>, Mixture<Huber, Huber>>* MX = new Mixture<Mixture<Huber, Huber>, Mixture<Huber, Huber>>(MX1, MX2, 0.5);


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
	cout << "Математическое ожидание: " << MX->M_Ksi() << endl;
	cout << "Дисперсия: " << MX->D_Ksi() << endl;
	cout << "Коэффиицент асимметрии: " << MX->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX->kurtosis() << endl << endl;

	Empirical* EM1 = new Empirical(MX, 1000, 0);
	Empirical* EM2 = new Empirical(MX, 1000, 0);

	cout << "Эмпирические характеристики:" << endl;
	cout << "Математическое ожидание: " << MX->M_Ksi() << endl;
	cout << "Дисперсия: " << MX->D_Ksi() << endl;
	cout << "Коэффициент асимметрии: " << MX->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX->kurtosis() << endl << endl;


	vector<pair<double, double>> s1 = MX->generate_pair(1000, MX->selection(1000));
	vector<pair<double, double>> s2 = EM1->generate_pair(1000, EM1->selection(1000));
	vector<pair<double, double>> s3 = EM2->generate_pair(1000, EM2->selection(1000));

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


void demo()
{
	cout << "1. Стандартная ситуация:" << endl << endl;

	Huber* HB1 = new Huber(1, 2, 0);
	Huber* HB2 = new Huber(1, 8, 0);

	Mixture<Huber, Huber>* MX1 = new Mixture<Huber, Huber>(HB1, HB2, 0.5);

	cout << "Параметр формы 2: " << MX1->get_component1()->get_v() << endl;
	cout << "Параметр масштаба 2: " << MX1->get_component1()->get_scale() << endl;
	cout << "Параметр сдвига 2: " << MX1->get_component1()->get_shift() << endl << endl;

	cout << "Параметр формы 2: " << HB2->get_v() << endl;
	cout << "Параметр масштаба 2: " << HB2->get_scale() << endl;
	cout << "Параметр сдвига 2: " << HB2->get_shift() << endl << endl;

	cout << "Теоретические характеристики:" << endl;
	cout << "Математическое ожидание: " << MX1->M_Ksi() << endl;
	cout << "Дисперсия: " << MX1->D_Ksi() << endl;
	cout << "Коэффиицент асимметрии: " << MX1->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX1->kurtosis() << endl << endl;

	Empirical* EM1 = new Empirical(HB1, 1000, 0);
	Empirical* EM2 = new Empirical(HB2, 1000, 0);

	cout << "Эмпирические характеристики:" << endl;
	cout << "Математическое ожидание: " << MX1->M_Ksi() << endl;
	cout << "Дисперсия: " << MX1->D_Ksi() << endl;
	cout << "Коэффициент асимметрии: " << MX1->asymmetry() << endl;
	cout << "Коэффициент эксцесса: " << MX1->kurtosis() << endl << endl;

	vector<pair<double, double>> s1 = MX1->generate_pair(1000, MX1->selection(1000));
	vector<pair<double, double>> s2 = EM1->generate_pair(1000, EM1->selection(1000));
	vector<pair<double, double>> s3 = EM2->generate_pair(1000, EM2->selection(1000));

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
}

int main(int argc, char** argv)
{
	setlocale(LC_ALL, "ru");

	demo();

	//start();

	//int result = Catch::Session().run(argc, argv);
	//return result;

	return 0;
}