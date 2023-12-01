#include "catch.hpp"
#include "function.h"
#include "function.cpp"
#include "Distribution.h"


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

    CHECK(HD->density(0) == Approx(0.34).epsilon(0.01));
    CHECK(HD->M_Ksi() == Approx(0).epsilon(0.01));
    CHECK(HD->D_Ksi() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Shift Scale Transformation")
{
    HuberD* HD = new HuberD();
    HD->set_scale(2);
    HD->set_shift(2);

    CHECK(HD->density(0) == Approx(0.103).epsilon(0.01));
    CHECK(HD->M_Ksi() == Approx(2).epsilon(0.01));
    CHECK(HD->D_Ksi() == Approx(2.24).epsilon(0.01));
    CHECK(HD->asymmetry() == Approx(0).epsilon(0.01));
    CHECK(HD->kurtosis() == Approx(2.37).epsilon(0.01));
}

TEST_CASE("Mixture Distribution")
{
    HuberD* HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(2);
    HD2->set_scale(2);
    HD1->set_shift(2);
    HD2->set_shift(2);
    Mixture<HuberD, HuberD>* MD = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);

    CHECK(MD->density(0) == Approx(0.103).epsilon(0.01));
    CHECK(MD->M_Ksi() == Approx(2).epsilon(0.01));
    CHECK(MD->D_Ksi() == Approx(2.24).epsilon(0.01));
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
    Mixture<HuberD, HuberD>* MD = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);

    CHECK(MD->M_Ksi() == Approx(1.5).epsilon(0.01));//MKsi_mixture
}

TEST_CASE("Mixture Distribution - Variance Test")
{
    HuberD* HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(1);
    HD2->set_scale(3);
    Mixture<HuberD, HuberD>* MD = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);

    CHECK(MD->Dksi_Mixture() == Approx(2.24).epsilon(0.01));//?
}

TEST_CASE("Late Binding Mechanism")
{
    HuberD* HD1 = new HuberD();
    HuberD* HD2 = new HuberD();
    HD1->set_scale(1);
    HD2->set_scale(3);

    Mixture<HuberD, HuberD>* MD1 = new Mixture<HuberD, HuberD>(HD1, HD2, 0.5);
    Mixture<HuberD, Mixture<HuberD, HuberD>>* MD2 = new Mixture<HuberD, Mixture<HuberD, HuberD>>(HD1, MD1, 0.5);

    CHECK(MD2->get_component1()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MD2->get_component1()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MD2->get_component1()->get_scale() == Approx(1).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component1()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component1()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component1()->get_scale() == Approx(1).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component2()->get_v() == Approx(1).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component2()->get_shift() == Approx(0).epsilon(0.01));
    CHECK(MD2->get_component2()->get_component2()->get_scale() == Approx(3).epsilon(0.01));
}

TEST_CASE("Empirical Distribution")
{
    HuberD* HD1 = new HuberD();
    Empirical* ED = new Empirical(HD1, 200, 1);

    CHECK(ED->get_n() == 200);
    CHECK(ED->get_k() == 8);
}