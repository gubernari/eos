/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2026 Nico Gubernari
 *
 * This file is part of the EOS project. EOS is free software;
 * you can redistribute it and/or modify it under the terms of the GNU General
 * Public License version 2, as published by the Free Software Foundation.
 *
 * EOS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 59 Temple
 * Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include <test/test.hh>
#include <eos/form-factors/parametric-g2026-impl.hh>

using namespace test;
using namespace eos;

class BToKG2026FormFactorsTest :
    public TestCase
{
    public:
        BToKG2026FormFactorsTest() :
            TestCase("b_to_k_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K::a^f+_0@G2026"]  =  0.01;
                p["B->K::a^f+_1@G2026"]  = -0.02;
                p["B->K::a^f0_1@G2026"]  =  0.05;
                p["B->K::a^fT_0@G2026"]  =  0.03;
                p["B->K::a^fT_1@G2026"]  = -0.04;

                p["mass::B_d@BSZ2015"]   =  5.279;
                p["mass::K_d@BSZ2015"]   =  0.494;

                p["mass::B_s,A^0[1]@G2026"] = 5.367;
                p["mass::B_s,V^0[1]@G2026"] = 5.711;
                p["mass::B_s,V^1[1]@G2026"] = 5.415;
                p["mass::B_s,A^1[1]@G2026"] = 5.829;

                // Optimized s0 = (mB + mK) * (sqrt(mB) - sqrt(mK))^2
                p["B->K::s0@G2026"]      =  14.682165;
                p["B->K::sV@G2026"]      =  30.272004;
                p["B->K::Q2@G2026"]      =  0.0;

                p["B->K::tchi_V0@G2026"]  =  1.42e-2;
                p["B->K::tchi_V1@G2026"]  =  7.30e-4;
                p["B->K::tchi_T1@G2026"]  =  4.55e-4;

                G2026FormFactors<BToK, PToP> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.16439553,  eps), // z(q2 =  0)
                    std::make_pair(  0.052958941, eps), // z(q2 = 11)
                    std::make_pair( -0.059753422, eps), // z(q2 = 18)
                    std::make_pair(  0.18728718,  eps), // z(q2 = -3)

                    std::make_pair(  1.0,           eps), // z_0(z = 0)
                    std::make_pair(  0.16439553,    eps), // z_1(z = 0)
                    std::make_pair(  0.027025890,   eps), // z_2(z = 0)
                    std::make_pair(  0.0044429356,  eps), // z_3(z = 0)
                    std::make_pair(  0.00073039876, eps), // z_4(z = 0)
                    std::make_pair(  0.00012007429, eps), // z_5(z = 0)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19))
                    std::make_pair( -0.080897132,    eps), // z_1(z = z(q2 = 19))
                    std::make_pair(  0.0065443460,   eps), // z_2(z = z(q2 = 19))
                    std::make_pair( -0.00052941882,  eps), // z_3(z = z(q2 = 19))
                    std::make_pair(  0.000042828465, eps), // z_4(z = z(q2 = 19))
                    std::make_pair( -0.0000034647,   eps), // z_5(z = z(q2 = 19))

                    std::make_pair(  1.0,           eps), // z_0(z = z(q2 = -3))
                    std::make_pair(  0.18728718,    eps), // z_1(z = z(q2 = -3))
                    std::make_pair(  0.035076488,   eps), // z_2(z = z(q2 = -3))
                    std::make_pair(  0.0065693766,  eps), // z_3(z = z(q2 = -3))
                    std::make_pair(  0.0012303600,  eps), // z_4(z = z(q2 = -3))
                    std::make_pair(  0.00023043066, eps), // z_5(z = z(q2 = -3))

                    std::make_pair(  0.09925017, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.09477547, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.08994813, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.06678005, eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.06447112, eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.06187060, eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.2435218,  eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.2269179,  eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.2097406,  eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  -0.0017120938,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.055553910, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.11501982,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 0.97754178,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.15219810,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.083367774, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0),-0.52688301,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.084602925, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.16545090,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 1.5162935,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V0(),       0.0025029313, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A0(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A1(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V1(),       0.0005,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_T1(),       0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_AT1(),      0,            eps);

                // Test everything for different s0 and sV
                p["B->K::s0@G2026"]      = -4.0;
                p["B->K::sV@G2026"]      =  33.327529;

                G2026FormFactors<BToK, PToP> ff2(p, Options{ });

                diagnostics = ff2.diagnostics();
                static const std::vector<std::pair<double, double>> reference2
                {
                    std::make_pair( -0.028329254,  eps), // z(q2 =  0)
                    std::make_pair( -0.12777540,   eps), // z(q2 = 11)
                    std::make_pair( -0.21891875,   eps), // z(q2 = 18)
                    std::make_pair( -0.0067887129, eps), // z(q2 = -3)

                    std::make_pair(  1.0,            eps), // z_0(z = 0)
                    std::make_pair( -0.028329254,    eps), // z_1(z = 0)
                    std::make_pair(  0.00080254664,  eps), // z_2(z = 0)
                    std::make_pair( -0.000022735548, eps), // z_3(z = 0)
                    std::make_pair(  6.4408111e-7,   eps), // z_4(z = 0)
                    std::make_pair( -1.8246337e-8,   eps), // z_5(z = 0)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19))
                    std::make_pair( -0.23491673,     eps), // z_1(z = z(q2 = 19))
                    std::make_pair(  0.055185872,    eps), // z_2(z = z(q2 = 19))
                    std::make_pair( -0.012964085,    eps), // z_3(z = z(q2 = 19))
                    std::make_pair(  0.0030454805,   eps), // z_4(z = z(q2 = 19))
                    std::make_pair( -0.00071543433,  eps), // z_5(z = z(q2 = 19))

                    std::make_pair(  1.0,            eps), // z_0
                    std::make_pair( -0.0067887129,   eps), // z_1
                    std::make_pair(  0.000046086623, eps), // z_2
                    std::make_pair( -3.1286885e-7,   eps), // z_3
                    std::make_pair(  2.1239768e-9,   eps), // z_4
                    std::make_pair( -1.4419069e-11,  eps), // z_5

                    std::make_pair(  0.09239401, eps), // phi_f_p(z = z(q2 = -2))
                    std::make_pair(  0.08860722, eps), // phi_f_p(z = z(q2 =  1))
                    std::make_pair(  0.08450451, eps), // phi_f_p(z = z(q2 =  4))
                    std::make_pair(  0.2153878,  eps), // phi_f_0(z = z(q2 = -2))
                    std::make_pair(  0.2172881,  eps), // phi_f_0(z = z(q2 =  1))
                    std::make_pair(  0.2190775,  eps), // phi_f_0(z = z(q2 =  4))
                    std::make_pair(  0.2375216,  eps), // phi_f_t(z = z(q2 = -2))
                    std::make_pair(  0.2227716,  eps), // phi_f_t(z = z(q2 =  1))
                    std::make_pair(  0.2074462,  eps), // phi_f_t(z = z(q2 =  4))

                    std::make_pair(  0.040519167,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(0.0),           ff2.f_p(0.0), eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.14861611, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.27516510, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 2.3634555,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.27005819, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.23639101, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.19320897, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.16757222, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.32309113, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 3.2127419,  eps);
            }
        }
} b_to_k_g2026_form_factors_test;

class BToKstarG2026FormFactorsTest :
    public TestCase
{
    public:
        BToKstarG2026FormFactorsTest() :
            TestCase("b_to_kstar_g2026_form_factors_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1e-5;

            {
                Parameters p = Parameters::Defaults();
                p["B->K^*::a^V_0@G2026"]     =  0.01;
                p["B->K^*::a^V_1@G2026"]     = -0.02;
                p["B->K^*::a^A0_0@G2026"]    =  0.03;
                p["B->K^*::a^A0_1@G2026"]    = -0.04;
                p["B->K^*::a^A1_1@G2026"]    =  0.05;
                p["B->K^*::a^A12_1@G2026"]   = -0.06;
                p["B->K^*::a^T1_0@G2026"]    =  0.07;
                p["B->K^*::a^T1_1@G2026"]    = -0.08;
                p["B->K^*::a^T2_1@G2026"]    =  0.09;
                p["B->K^*::a^T23_1@G2026"]   = -0.10;

                p["mass::B_d@BSZ2015"]     =  5.279;
                p["mass::K_d^*@BSZ2015"]   =  0.892;

                p["mass::B_s,A^0[1]@G2026"] = 5.367;
                p["mass::B_s,V^0[1]@G2026"] = 5.711;
                p["mass::B_s,V^1[1]@G2026"] = 5.415;
                p["mass::B_s,A^1[1]@G2026"] = 5.829;

                // Optimized s0 = (mB + mK) * (sqrt(mB) - sqrt(mK))^2
                p["B->K^*::s0@G2026"]      =  11.299192;
                p["B->K^*::sV@G2026"]      =  30.272004;
                p["B->K^*::sA@G2026"]      =  31.775769;
                p["B->K^*::Q2@G2026"]      =  0.0;

                p["B->K^*::tchi_A0@G2026"]  =  1.52e-2;
                p["B->K^*::tchi_V1@G2026"]  =  6.97e-4;
                p["B->K^*::tchi_A1@G2026"]  =  6.54e-4;
                p["B->K^*::tchi_T1@G2026"]  =  4.83e-4;
                p["B->K^*::tchi_AT1@G2026"] =  4.12e-4;

                G2026FormFactors<BToKstar, PToV> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.11627577,   eps), // z(q2 =  0, sV)
                    std::make_pair(  0.0039115965, eps), // z(q2 = 11, sV)
                    std::make_pair( -0.098769544,  eps), // z(q2 = 18, sA)
                    std::make_pair(  0.13164137,   eps), // z(q2 = -3, sA)

                    std::make_pair(  1.0,            eps), // z_0(z = 0, sV)
                    std::make_pair(  0.11627577,     eps), // z_1(z = 0, sV)
                    std::make_pair(  0.013520054,    eps), // z_2(z = 0, sV)
                    std::make_pair(  0.0015720547,   eps), // z_3(z = 0, sV)
                    std::make_pair(  0.00018279187,  eps), // z_4(z = 0, sV)
                    std::make_pair(  0.000021254265, eps), // z_5(z = 0, sV)

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = 19, sV))
                    std::make_pair( -0.12944094,     eps), // z_1(z = z(q2 = 19, sV))
                    std::make_pair(  0.016754956,    eps), // z_2(z = z(q2 = 19, sV))
                    std::make_pair( -0.0021687773,   eps), // z_3(z = z(q2 = 19, sV))
                    std::make_pair(  0.00028072856,  eps), // z_4(z = z(q2 = 19, sV))
                    std::make_pair( -0.000036337769, eps), // z_5(z = z(q2 = 19, sV))

                    std::make_pair(  1.0,            eps), // z_0(z = z(q2 = -3, sA))
                    std::make_pair(  0.13164137,     eps), // z_1(z = z(q2 = -3, sA))
                    std::make_pair(  0.017329451,    eps), // z_2(z = z(q2 = -3, sA))
                    std::make_pair(  0.0022812726,   eps), // z_3(z = z(q2 = -3, sA))
                    std::make_pair(  0.00030030986,  eps), // z_4(z = z(q2 = -3, sA))
                    std::make_pair(  0.000039533201, eps), // z_5(z = z(q2 = -3, sA))

                    std::make_pair(  0.3152006,  eps), // phi_v(z = z(q2 = -2))
                    std::make_pair(  0.2974732,  eps), // phi_v(z = z(q2 =  1))
                    std::make_pair(  0.2790677,  eps), // phi_v(z = z(q2 =  4))

                    std::make_pair(  0.5025545,  eps), // phi_a_0(z = z(q2 = -2))
                    std::make_pair(  0.4746616,  eps), // phi_a_0(z = z(q2 =  1))
                    std::make_pair(  0.4457290,  eps), // phi_a_0(z = z(q2 =  4))
                    std::make_pair(  0.05659679, eps), // phi_a_1(z = z(q2 = -2))
                    std::make_pair(  0.05492700, eps), // phi_a_1(z = z(q2 =  1))
                    std::make_pair(  0.05304411, eps), // phi_a_1(z = z(q2 =  4))
                    std::make_pair(  0.02133889, eps), // phi_a_12(z = z(q2 = -2))
                    std::make_pair(  0.02119834, eps), // phi_a_12(z = z(q2 =  1))
                    std::make_pair(  0.02099215, eps), // phi_a_12(z = z(q2 =  4))

                    std::make_pair(  0.2089456, eps), // phi_t_1(z = z(q2 = -2))
                    std::make_pair(  0.2020819, eps), // phi_t_1(z = z(q2 =  1))
                    std::make_pair(  0.1946578, eps), // phi_t_1(z = z(q2 =  4))

                    std::make_pair(  0.02732397, eps), // phi_t_2(z = z(q2 = -2))
                    std::make_pair(  0.02714399, eps), // phi_t_2(z = z(q2 =  1))
                    std::make_pair(  0.02687998, eps), // phi_t_2(z = z(q2 =  4))
                    std::make_pair(  0.02493920, eps), // phi_t_23(z = z(q2 = -2))
                    std::make_pair(  0.02420341, eps), // phi_t_23(z = z(q2 =  1))
                    std::make_pair(  0.02337372, eps), // phi_t_23(z = z(q2 =  4))

                    std::make_pair(  0.00808554091, eps), // a_A12_0
                    std::make_pair(  0.05171365,    eps), // a_A1_0
                    std::make_pair(  0.001710846,   eps), // a_T2_0
                    std::make_pair(  -0.02595125,   eps), // a_T23_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                const double mB = p["mass::B_d@BSZ2015"];
                const double mV = p["mass::K_d^*@BSZ2015"];
                const double tm = (mB - mV) * (mB - mV);

                const double factora12a0 = (mB * mB - mV * mV) / 8.0 / mB / mV;
                const double factora12a1 = (mB + mV) * (mB * mB - mV * mV - tm) / 16.0 / mB / mV / mV;
                const double factort23t2 = (mB + mV) * (mB * mB + 3.0 * mV * mV - tm) / 8.0 / mB / mV / mV;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(tm) , factora12a1 * ff.a_1(tm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(tm) , factort23t2 * ff.t_2(tm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.019856051, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.041885252, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.33524673,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.059631278, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.11372013,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.0100095,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.99821129,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  1.0420689,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  1.3310723,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.19263127,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.14210536,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  1.6126605,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.30288816,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  0.46431383,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  2.0639932,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  0.72944768,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  0.34650161,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0), -1.1708205,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0), -1.6950877,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0), -1.4553959,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  0.080068039, eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V0(),  0,          eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A0(),  0.0025,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A1(),  0.00883968, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V1(),  0.0005,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_T1(),  0.0113,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_AT1(), 0.0187764,  eps);


                // Test everything for different s0 and sV
                p["B->K^*::s0@G2026"] = -4.0;
                p["B->K^*::sV@G2026"] =  38.081241;
                p["B->K^*::sA@G2026"] =  38.081241;

                // Test end-point relations

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(tm) , factora12a1 * ff.a_1(tm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(tm) , factort23t2 * ff.t_2(tm) , eps);

                // Test against my Mathematica implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   (-15.0),  0.060354646, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  3.0),  0.12451300,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  1.2368289,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (-15.0),  0.12014397,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  3.0),  0.23335484,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  2.3503654,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (-15.0),  0.69051914,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  3.0),  0.69204300,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.93515203,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(-15.0), -0.021451283, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  3.0),  0.18777326,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.92333822,  eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (-15.0),  0.81003891,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  3.0),  1.3711616,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  9.6563017,   eps);

                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (-15.0),  1.4521739,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  3.0),  1.1807589,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.85018484,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(-15.0),  0.67859760,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  3.0),  1.0230491,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  2.8515852,   eps);


                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V0(),  0,          eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A0(),  0.0025,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A1(),  0.00969611, eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V1(),  0.0005,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_T1(),  0.0113,     eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_AT1(), 0.0215922,  eps);
            }
        }
} b_to_kstar_g2026_form_factors_test;
