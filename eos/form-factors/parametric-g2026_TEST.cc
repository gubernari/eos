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

                p["mass::B_s,A^0@G2026"] =  5.367;
                p["mass::B_s,V^0*@G2026"] = 5.711;
                p["mass::B_s,V^1@G2026"] =  5.415;
                p["mass::B_s,A^1*@G2026"] = 5.829;

                // Optimized s0 = (mB + mK) * (sqrt(mB) - sqrt(mK))^2
                p["B->K::s0@G2026"]      =  14.682165;
                p["B->K::sV@G2026"]      =  30.272004;
                p["B->K::Q2@G2026"]      =  0.0;

                p["B->K::tchi_V0@G2026"]  =  1.42e-2;
                p["B->K::tchi_V1@G2026"]  =  7.30e-4;
                p["B->K::tchi_T1@G2026"]  =  4.55e-4;

                G2026FormFactors<BToK, PToP> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                /*static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.16439553,  eps), // z(q2 =  0)
                    std::make_pair(  0.052958941, eps), // z(q2 = 11)
                    std::make_pair( -0.059753422, eps), // z(q2 = 18)
                    std::make_pair(  0.18728718,  eps), // z(q2 = -3)

                    std::make_pair(  1.0,           eps), // z_0(z = 0.0)
                    std::make_pair(  0.16439553,    eps), // z_1(z = 0.0)
                    std::make_pair(  0.027025890,   eps), // z_2(z = 0.0)
                    std::make_pair(  0.0044429356,  eps), // z_3(z = 0.0)
                    std::make_pair(  0.00073039876, eps), // z_4(z = 0.0)
                    std::make_pair(  0.00012007429, eps), // z_5(z = 0.0)

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
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);*/

                // Test end-point relations
                //TEST_CHECK_NEARLY_EQUAL( ff.f_0(0.0),           ff.f_p(0.0), eps);

                // Test against Nico's implementation
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.055553910, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.11501982,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 0.97754178,  eps);
                //TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.15219810,  eps);
                //TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.083367774, eps);
                //TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0),-0.52688301,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.084602925, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.16545090,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 1.5162935,   eps);

                //TEST_CHECK_NEARLY_EQUAL( ff.saturation_V0(),       0.0025029313, eps);
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

                    std::make_pair(  1.0,            eps), // z_0(z = 0.0)
                    std::make_pair( -0.028329254,    eps), // z_1(z = 0.0)
                    std::make_pair(  0.00080254664,  eps), // z_2(z = 0.0)
                    std::make_pair( -0.000022735548, eps), // z_3(z = 0.0)
                    std::make_pair(  6.4408111e-7,   eps), // z_4(z = 0.0)
                    std::make_pair( -1.8246337e-8,   eps), // z_5(z = 0.0)

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

                    std::make_pair(  -0.0017120938,  eps), // a_f0_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference2);

                // Test end-point relations
                TEST_CHECK_NEARLY_EQUAL( ff2.f_0(0.0),           ff2.f_p(0.0), eps);

                // Test against Nico's implementation
                //TEST_CHECK_RELATIVE_ERROR( ff.f_p(-15.0), 0.13869773, eps);
                //TEST_CHECK_RELATIVE_ERROR( ff.f_p(  3.0), 0.26133659, eps);
                //TEST_CHECK_RELATIVE_ERROR( ff.f_p( 25.0), 2.3173841,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(-15.0), 0.25648978, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0(  3.0), 0.22355027, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_0( 25.0), 0.17554911, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(-15.0), 0.16033129, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t(  3.0), 0.31202406, eps);
                TEST_CHECK_RELATIVE_ERROR( ff.f_t( 25.0), 3.1713797,  eps);
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

                p["mass::B_d@BSZ2015"]        =  5.279;
                p["mass::K_d^*@BSZ2015"]      =  0.896;
                p["mass::B_s@BSZ2015"]        =  5.367;
                p["mass::B_s^*@BSZ2015"]      =  5.416;
                p["mass::B_s,0@BSZ2015"]      =  5.711;
                p["mass::B_s,1@BSZ2015"]      =  5.750;

                // Optimized s0 = (mB + mK*) * (sqrt(mB) - sqrt(mK*))^2
                p["B->K^*::s0@G2026"]       =  11.271194912;
                p["B->K^*::sV@G2026"]       =  30.261001;
                p["B->K^*::sA@G2026"]       =  31.764496;
                p["B->K^*::Q2@G2026"]       =  0.0;

                G2026FormFactors<BToKstar, PToV> ff(p, Options{ });

                Diagnostics diagnostics = ff.diagnostics();
                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(  0.109126,  eps), // z_a(q2 =  0)
                    std::make_pair(  0.115965,  eps), // z_v(q2 =  0)
                    std::make_pair(  0.015044,  eps), // z_a(q2 = 10)
                    std::make_pair(  0.016197,  eps), // z_v(q2 = 10)
                    std::make_pair(  0.500293,  eps), // p_0(z = 0.0)
                    std::make_pair( -0.256101,  eps), // p_1(z = 0.0)
                    std::make_pair(  0.324605,  eps), // p_2(z = 0.0)
                    std::make_pair( -0.395358,  eps), // p_3(z = 0.0)
                    std::make_pair(  0.474303,  eps), // p_4(z = 0.0)
                    std::make_pair( -0.565749,  eps), // p_5(z = 0.0)
                    std::make_pair(  0.500293,  eps), // p_0(z = z(q2 = 10))
                    std::make_pair( -0.246998,  eps), // p_1(z = z(q2 = 10))
                    std::make_pair(  0.317589,  eps), // p_2(z = z(q2 = 10))
                    std::make_pair( -0.385009,  eps), // p_3(z = z(q2 = 10))
                    std::make_pair(  0.459807,  eps), // p_4(z = z(q2 = 10))
                    std::make_pair( -0.545962,  eps), // p_5(z = z(q2 = 10))

                    std::make_pair(  0.127438, eps), // phi_v(z = z(q2 = -2.0))
                    std::make_pair(  0.120283, eps), // phi_v(z = z(q2 =  1.0))
                    std::make_pair(  0.112854, eps), // phi_v(z = z(q2 =  4.0))
                    std::make_pair(  0.197626, eps), // phi_a_0(z = z(q2 = -2.0))
                    std::make_pair(  0.186676, eps), // phi_a_0(z = z(q2 =  1.0))
                    std::make_pair(  0.175318, eps), // phi_a_0(z = z(q2 =  4.0))
                    std::make_pair(  0.083246, eps), // phi_a_1(z = z(q2 = -2.0))
                    std::make_pair(  0.084125, eps), // phi_a_1(z = z(q2 =  1.0))
                    std::make_pair(  0.084999, eps), // phi_a_1(z = z(q2 =  4.0))
                    std::make_pair(  0.031512, eps), // phi_a_12(z = z(q2 = -2.0))
                    std::make_pair(  0.032597, eps), // phi_a_12(z = z(q2 =  1.0))
                    std::make_pair(  0.033773, eps), // phi_a_12(z = z(q2 =  4.0))
                    std::make_pair(  0.086038, eps), // phi_t_1(z = z(q2 = -2.0))
                    std::make_pair(  0.083221, eps), // phi_t_1(z = z(q2 =  1.0))
                    std::make_pair(  0.080174, eps), // phi_t_1(z = z(q2 =  4.0))
                    std::make_pair(  0.039178, eps), // phi_t_2(z = z(q2 = -2.0))
                    std::make_pair(  0.040526, eps), // phi_t_2(z = z(q2 =  1.0))
                    std::make_pair(  0.041989, eps), // phi_t_2(z = z(q2 =  4.0))
                    std::make_pair(  0.035899, eps), // phi_t_23(z = z(q2 = -2.0))
                    std::make_pair(  0.036278, eps), // phi_t_23(z = z(q2 =  1.0))
                    std::make_pair(  0.036654, eps), // phi_t_23(z = z(q2 =  4.0))

                    std::make_pair(  0.10207 , eps), // a_A1_0
                    std::make_pair( -0.009349, eps), // a_A12_0
                    std::make_pair(  0.098888, eps), // a_T2_0
                    std::make_pair(  0.013503, eps), // a_T23_0
                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                // Test end-point relations
                const double tm = (BToKstar::m_B - BToKstar::m_V) * (BToKstar::m_B - BToKstar::m_V);

                const double factora12a0 = (BToKstar::m_B * BToKstar::m_B - BToKstar::m_V * BToKstar::m_V) / 8.0 / BToKstar::m_B / BToKstar::m_V;
                const double factora12a1 = (BToKstar::m_B + BToKstar::m_V) * (BToKstar::m_B * BToKstar::m_B - BToKstar::m_V * BToKstar::m_V - tm)
                                           / 16.0 / BToKstar::m_B / BToKstar::m_V / BToKstar::m_V;
                const double factort23t2 = (BToKstar::m_B + BToKstar::m_V) * (BToKstar::m_B * BToKstar::m_B + 3.0 * BToKstar::m_V * BToKstar::m_V - tm)
                                           / 8.0 / BToKstar::m_B / BToKstar::m_V / BToKstar::m_V;

                TEST_CHECK_NEARLY_EQUAL( ff.a_12(0.0), factora12a0 * ff.a_0(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_1(0.0) ,               ff.t_2(0.0), eps);
                TEST_CHECK_NEARLY_EQUAL( ff.a_12(tm) , factora12a1 * ff.a_1(tm) , eps);
                TEST_CHECK_NEARLY_EQUAL( ff.t_23(tm) , factort23t2 * ff.t_2(tm) , eps);

                // Test against Nico's implementation
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( -1.0),  0.098866,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  1.0),  0.10623,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   (  4.0),  0.119471,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.v   ( 25.0),  0.688017,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( -1.0),  0.197014,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  1.0),  0.212504,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 (  4.0),  0.2407,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_0 ( 25.0),  1.80208,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( -1.0),  0.502886,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  1.0),  0.494491,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 (  4.0),  0.48148,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_1 ( 25.0),  0.361849,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( -1.0),  0.140285,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  1.0),  0.152292,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12(  4.0),  0.170921,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.a_12( 25.0),  0.362945,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( -1.0),  0.830824,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  1.0),  0.87265,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 (  4.0),  0.946545,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_1 ( 25.0),  3.62504,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( -1.0),  0.869941,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  1.0),  0.8321,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 (  4.0),  0.774239,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_2 ( 25.0),  0.2791,    eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( -1.0),  0.59847,   eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  1.0),  0.617008,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23(  4.0),  0.647415,  eps);
                TEST_CHECK_RELATIVE_ERROR( ff.t_23( 25.0),  1.11061,   eps);

                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V0(),       0,            eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A0(),       0.0025,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_A1(),       0.0005,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_V1(),       0.0166056,    eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_T1(),       0.0113,       eps);
                TEST_CHECK_NEARLY_EQUAL( ff.saturation_AT1(),      0.0280612,    eps);
            }
        }
} b_to_kstar_g2026_form_factors_test;
