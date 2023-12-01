#include "catch.hpp"
#include "function.h"


TEST_CASE("Basic Methods")
{
    HuberD* HD = new HuberD();

    CHECK(HD->get_v() == Approx(1).epsilon(0.01));
    CHECK(HD->get_shift() == Approx(0).epsilon(0.01));
    CHECK(HD->get_scale() == Approx(1).epsilon(0.01));
}

TEST_CASE("Standard Distribution")
{
    HuberD* HD = new HuberD();

    CHECK(HD->Huber(0) == Approx(0.34).epsilon(0.01));
    CHECK(HD->Mksi_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->Dksi_huber() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis_huber() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Shift Scale Transformation")
{
    HuberD* HD = new HuberD();
    HD->set_scale(2);
    HD->set_shift(2);

    CHECK(HD->Huber(0) == Approx(0.103).epsilon(0.01));
    CHECK(HD->Mksi_huber() == Approx(2).epsilon(0.01));
    CHECK(HD->Dksi_huber() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry_huber() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis_huber() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Mixture Distribution")
{
    HuberD * HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(2);
    HD2->set_scale(2);
    HD1->set_shift(2);
    HD2->set_shift(2);
    Mixture<HuberD, HuberD>* MD = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);

    CHECK(MD->H_Mixture(0) == Approx(0.103).epsilon(0.01));
    CHECK(MD->expected_value() == Approx(2).epsilon(0.01));
    CHECK(MD->variance() == Approx(2.24).epsilon(0.01));
    CHECK(MD->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(MD->kurtosis() == Approx(1.77).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Expected Value Test")
{
    HuberD* HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(2);
    HD2->set_scale(2);
    HD1->set_shift(1);
    HD2->set_shift(2);
    Mixture<Primary, Primary>* MX = new Mixture<Primary, Primary>(HB1, HB2, 0.5);

    CHECK(MX->expected_value() == Approx(1.5).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Variance Test")
{
    HuberD* HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(1);
    HD2->set_scale(3);
    Mixture<Primary, Primary>* MX = new Mixture<Primary, Primary>(HB1, HB2, 0.5);

    CHECK(MX->variance() == Approx(2.24).epsilon(0.01));
}

TEST_CASE("Late Binding Mechanism")
{
    Primary* HB1 = new Primary();
    Primary* HB2 = new Primary();
    HB1->set_scale(1);
    HB2->set_scale(3);

    Mixture<Primary, Primary>* MX1 = new Mixture<Primary, Primary>(HB1, HB2, 0.5);
    Mixture<Primary, Mixture<Primary, Primary>>* MX2 = new Mixture<Primary, Mixture<Primary, Primary>>(HB1, MX1, 0.5);

    CHECK(MX2->get_component1()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MX2->get_component1()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MX2->get_component1()->get_scale() == Approx(1).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component1()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component1()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component1()->get_scale() == Approx(1).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component2()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component2()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MX2->get_component2()->get_component2()->get_scale() == Approx(3).epsilon(0.01));
}

TEST_CASE("Empirical Distribution")
{
    Primary* HB1 = new Primary();
    Empirical* EM = new Empirical(HB1, 200, 1);

    CHECK(EM->get_n() == 200);
    CHECK(EM->get_k() == 8);
}