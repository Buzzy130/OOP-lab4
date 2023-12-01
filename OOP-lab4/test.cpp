#include "catch.hpp"
#include "function.h"
#include "Mixture.cpp"


TEST_CASE("Basic Methods")
{
	Huber* HB = new Huber();

    CHECK(HB->get_v() == Approx(1).epsilon(0.01));
    CHECK(HB->get_shift() == Approx(0).epsilon(0.01));
    CHECK(HB->get_scale() == Approx(1).epsilon(0.01));
}

TEST_CASE("Standard Distribution")
{
    Huber* HB = new Huber();

    CHECK(HB->density(0) == Approx(0.34).epsilon(0.01));
    CHECK(HB->M_Ksi() == Approx(0).epsilon(0.01));
    CHECK(HB->D_Ksi() == Approx(2.24).epsilon(0.01));
    CHECK(HB->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HB->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Shift Scale Transformation")
{
    Huber* HB = new Huber();
    HB->set_scale(2);
    HB->set_shift(2);

    CHECK(HB->density(0) == Approx(0.103).epsilon(0.01));
    CHECK(HB->M_Ksi() == Approx(2).epsilon(0.01));
    CHECK(HB->D_Ksi() == Approx(2.24).epsilon(0.01));
    CHECK(HB->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HB->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Mixture Distribution")
{
    Huber* HB1 = new Huber();
    Huber* HB2 = new Huber();
    HB1->set_scale(2);
    HB2->set_scale(2);
    HB1->set_shift(2);
    HB2->set_shift(2);
    Mixture<Huber, Huber>* MX = new Mixture<Huber, Huber>(HB1, HB2, 0.5);

    CHECK(MX->density(0) == Approx(0.103).epsilon(0.01));
    CHECK(MX->M_Ksi() == Approx(2).epsilon(0.01));
    CHECK(MX->D_Ksi() == Approx(2.24).epsilon(0.01));
    CHECK(MX->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(MX->kurtosis() == Approx(1.77).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Expected Value Test")
{
    Huber* HB1 = new Huber();
    Huber* HB2 = new Huber();
    HB1->set_scale(2);
    HB2->set_scale(2);
    HB1->set_shift(1);
    HB2->set_shift(2);
    Mixture<Huber, Huber>* MX = new Mixture<Huber, Huber>(HB1, HB2, 0.5);

    CHECK(MX->M_Ksi() == Approx(1.5).epsilon(0.01));
}

TEST_CASE("Mixture Distribution - Variance Test")
{
    Huber* HB1 = new Huber();
    Huber* HB2 = new Huber();
    HB1->set_scale(1);
    HB2->set_scale(3);
    Mixture<Huber, Huber>* MX = new Mixture<Huber, Huber>(HB1, HB2, 0.5);

    CHECK(MX->D_Ksi() == Approx(2.24).epsilon(0.01));
}

TEST_CASE("Late Binding Mechanism")
{
    Huber* HB1 = new Huber();
    Huber* HB2 = new Huber();
    HB1->set_scale(1);
    HB2->set_scale(3);

    Mixture<Huber, Huber>* MX1 = new Mixture<Huber, Huber>(HB1, HB2, 0.5);
    Mixture<Huber, Mixture<Huber, Huber>>* MX2 = new Mixture<Huber, Mixture<Huber, Huber>>(HB1, MX1, 0.5);

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
    Huber* HB1 = new Huber();
    Empirical* EM = new Empirical(HB1, 200, 1);

    CHECK(EM->get_n() == 200);
    CHECK(EM->get_k() == 8);
}