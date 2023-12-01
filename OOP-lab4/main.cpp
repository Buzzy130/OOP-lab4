#define CATCH_CONFIG_RUNNER
#include "catch.hpp"
#include "function.h"
#include "Distribution.h"
#include "function.cpp"
#define _USE_MATH_DEFINES
#include <cmath>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <numeric>
#include <algorithm>
#include <string>
using namespace std;

void demo()
{
	cout << "7. ���� �������� ����������:" << endl << endl;

	HuberD* HD1 = new HuberD(1, 1, 2);
	HuberD* HD2 = new HuberD(1, 2, 6);
	HuberD* HD3 = new HuberD(1, 2, -6);
	HuberD* HD4 = new HuberD(1, 1, -2);

	Mixture<HuberD, HuberD>* MD1 = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);
	Mixture<HuberD, HuberD>* MD2 = new Mixture<HuberD, HuberD>(HD3, HD4, 0.5);
	Mixture<Mixture<HuberD, HuberD>, Mixture<HuberD, HuberD>>* MD = new Mixture<Mixture<HuberD, HuberD>, Mixture<HuberD, HuberD>>(MD1, MD2, 0.5);


	cout << "�������� ����� 1: " << MD->get_component1()->get_component1()->get_v() << endl;
	cout << "�������� �������� 1: " << MD->get_component1()->get_component1()->get_scale() << endl;
	cout << "�������� ������ 1: " << MD->get_component1()->get_component1()->get_shift() << endl << endl;

	cout << "�������� ����� 2: " << MD->get_component1()->get_component2()->get_v() << endl;
	cout << "�������� �������� 2: " << MD->get_component1()->get_component2()->get_scale() << endl;
	cout << "�������� ������ 2: " << MD->get_component1()->get_component2()->get_shift() << endl << endl;

	cout << "�������� ����� 3: " << MD->get_component2()->get_component1()->get_v() << endl;
	cout << "�������� �������� 3: " << MD->get_component2()->get_component1()->get_scale() << endl;
	cout << "�������� ������ 3: " << MD->get_component2()->get_component1()->get_shift() << endl << endl;

	cout << "�������� ����� 4: " << MD->get_component2()->get_component2()->get_v() << endl;
	cout << "�������� �������� 4: " << MD->get_component2()->get_component2()->get_scale() << endl;
	cout << "�������� ������ 4: " << MD->get_component2()->get_component2()->get_shift() << endl << endl;

	cout << "������������� ��������������:" << endl;
	cout << "�������������� ��������: " << MD->M_Ksi() << endl;//?
	cout << "���������: " << MD->D_Ksi() << endl;//?
	cout << "����������� ����������: " << MD->asymmetry() << endl;
	cout << "����������� ��������: " << MD->kurtosis() << endl << endl;

	Empirical* EM1 = new Empirical(MD, 10000, 0);
	Empirical* EM2 = new Empirical(MD, 10000, 0);

	cout << "������������ ��������������:" << endl;
	cout << "�������������� ��������: " << MD->M_Ksi() << endl;//?
	cout << "���������: " << MD->D_Ksi() << endl;//?
	cout << "����������� ����������: " << MD->asymmetry() << endl;
	cout << "����������� ��������: " << MD->kurtosis() << endl << endl;


	vector<pair<double, double>> s1 = MD->generate_pair(10000, MD->selection(10000));
	vector<pair<double, double>> s2 = EM1->generate_pair(10000, EM1->selection(10000));
	vector<pair<double, double>> s3 = EM2->generate_pair(10000, EM2->selection(10000));

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

	//demo();

	//start();

	int result = Catch::Session().run(argc, argv);
	return result;

	return 0;
}