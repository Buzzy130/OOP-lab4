#include "catch.hpp"
#include "function.h"
#include "function.cpp"

TEST_CASE("Basic Methods")
{
    Primary* HB = new Primary();

    CHECK(HB->get_v() == Approx(1).epsilon(0.01));
    CHECK(HB->get_shift() == Approx(0).epsilon(0.01));
    CHECK(HB->get_scale() == Approx(1).epsilon(0.01));
}

TEST_CASE("Standard Distribution")
{
    Primary* HB = new Primary();

    CHECK(HB->f(0) == Approx(0.34).epsilon(0.01));
    CHECK(HB->expected_value() == Approx(0).epsilon(0.01));
    CHECK(HB->variance() == Approx(2.24).epsilon(0.01));
    CHECK(HB->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HB->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Shift Scale Transformation")
{
    Primary* HB = new Primary();
    HB->set_scale(2);
    HB->set_shift(2);

    CHECK(HB->f(0) == Approx(0.103).epsilon(0.01));
    CHECK(HB->expected_value() == Approx(2).epsilon(0.01));
    CHECK(HB->variance() == Approx(2.24).epsilon(0.01));
    CHECK(HB->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HB->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Mixture Distribution")
{
    Primary* HB1 = new Primary();
    Primary* HB2 = new Primary();
    HB1->set_scale(2);
    HB2->set_scale(2);
    HB1->set_shift(2);
    HB2->set_shift(2);
    Mixture<Primary, Primary>* MX = new Mixture<Primary, Primary>(HB1, HB2, 0.5);

    CHECK(MX->f(0) == Approx(0.103).epsilon(0.01));
    CHECK(MX->expected_value() == Approx(2).epsilon(0.01));
    CHECK(MX->variance() == Approx(2.24).epsilon(0.01));
    CHECK(MX->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(MX->kurtosis() == Approx(1.77).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Expected Value Test")
{
    Primary* HB1 = new Primary();
    Primary* HB2 = new Primary();
    HB1->set_scale(2);
    HB2->set_scale(2);
    HB1->set_shift(1);
    HB2->set_shift(2);
    Mixture<Primary, Primary>* MX = new Mixture<Primary, Primary>(HB1, HB2, 0.5);

    CHECK(MX->expected_value() == Approx(1.5).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Variance Test")
{
    Primary* HB1 = new Primary();
    Primary* HB2 = new Primary();
    HB1->set_scale(1);
    HB2->set_scale(3);
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