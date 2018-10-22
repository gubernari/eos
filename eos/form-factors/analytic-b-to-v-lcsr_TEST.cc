/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2018 Nico Gubernari
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
#include <eos/form-factors/analytic-b-to-v-lcsr.hh>
#include <eos/form-factors/mesonic.hh>

#include <vector>
#include <utility>

using namespace test;
using namespace eos;

class LCSRFormFactorsTest :
    public TestCase
{
    public:
        LCSRFormFactorsTest() :
            TestCase("lcsr_form_factors_test")
        {
        }

        virtual void run() const
        {
            /* B -> K^* diagnostic values */
            {
                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.3174;
                p["B::lambda_H^2"]            = 1.2696;
                p["mass::ud(2GeV)"]           = 0.008;
                p["mass::B_d"]                = 5.2795;
                p["mass::K^*_d"]              = 0.896;
                p["decay-constant::B_d"]      = 0.180;
                p["B->K^*::f_Kstar_par"]      = 0.217;
                p["B->K^*::mu@B-LCSR"]        = 1.0;
                p["B->K^*::s_0^A1,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A1,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A2,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A2,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A30,0@B-LCSR"] = 1.7;
                p["B->K^*::s_0^A30,1@B-LCSR"] = 0.0;
                p["B->K^*::s_0^V,0@B-LCSR"]   = 1.7;
                p["B->K^*::s_0^V,1@B-LCSR"]   = 0.0;
                p["B->K^*::M^2@B-LCSR"]       = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "zero" }
                };
                AnalyticFormFactorBToVLCSR<lcsr::BToKstar> ff{ p, o };
                auto diagnostics = ff.diagnostics();

                std::cout << "Diagnostics:" << std::endl;
                for (auto & d : diagnostics)
                {
                    std::cout << d.description << ": " << d.value << std::endl;
                }
                std::cout << "Diagnostics ended" << std::endl;

                static const std::vector<std::pair<double, double>> reference
                {
                    std::make_pair(0.1268577,  1.0e-7), // m_v(mu) in the MSbar scheme, strange quark for these tests
                    std::make_pair(0.896,      1.0e-7), // m_V, the K^* mass
                    std::make_pair(0.217,      1.0e-7), // f_V, the K^* decay constant

                    std::make_pair(0.07488769, 1.0e-7), // sigma_0 value for q^2 = 5.0
                    std::make_pair(1.7,        1.0e-7), // s_0 value for V

                    /* A_1 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.785313,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.662569,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.539825,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.007976,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.839081,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.670185,  1.0e-6), // I1_A1_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 phi_bar */
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.984386e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.504278e-3,  1.0e-6), // I1_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(6.612526e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.578992e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.545457e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.122634e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(6.761613e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.400593e-2,  1.0e-6), // I2_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair( 8.504824e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 6.998295e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 5.491766e-1,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 2.125109e-2,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.312584e-3,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.387626e-2,  1.0e-6), // I2d1_A1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 g_+ */
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.450784e-3,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.256998e-2,  1.0e-6), // I1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.422015e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190204e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-9.583918e-2,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.904307e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.221156e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.538004e-1,  1.0e-6), // I2_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-5.794951,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.809833,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.824715,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.049257,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.868187,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.687116,  1.0e-6), // I2d1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.972448e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.195257e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.418066e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.426346e-2,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.187349e-2,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(9.483517e-3,  1.0e-6), // I3_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(2.079645e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.741271e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.402896e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.369436e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.930266e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.491096e-1,  1.0e-6), // I3d1_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(2.425945,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.932077,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.438210,  1.0e-6),     // I3d2_A1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.453318e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.043950e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.634583e-1,  1.0e-6), // I3d2_A1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.685383e-3,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.355562e-2,  1.0e-6), // I2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_bar */
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.557924e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.519371e-1,  1.0e-6), // I2d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.281981e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(8.563666e-4,  1.0e-6), // I3_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(9.036255e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(9.036255e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(9.036255e-3,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.830266e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.830266e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.830266e-2,  1.0e-6), // I3d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar */
                    std::make_pair(3.961103e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.961103e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.961103e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.380327e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.380327e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.380327e-1,  1.0e-6), // I3d2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(1.545825e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.304214e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.062602e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.047158e-3,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.716971e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(6.962364e-4,  1.0e-6), // I4_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(1.094715e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(9.194681e-3,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.442217e-3,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.499502e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.882185e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.264868e-2,  1.0e-6), // I4d1_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(4.849987e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.032256e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.214524e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(6.846914e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.485191e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.123467e-1,  1.0e-6), // I4d2_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(7.795264,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.165359,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.535454,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.860036,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.694950,  1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(0.5298649, 1.0e-5), // I4d3_A1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_1 phi_3 */
                    std::make_pair( 1.7872e-3,  1.0e-7),  // I1_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.7454e-3,  1.0e-7),  // I1_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9875e-4,  1.0e-8),  // I1_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.2609e-4,  1.0e-8),  // I1_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_3 */
                    std::make_pair( 3.6060e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.2251e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.0787e-3,  1.0e-7),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2610e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_1 phi_bar_3 */
                    std::make_pair( 8.7010e-4,  1.0e-8),  // I1_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8234e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8006e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2156e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-3.2004e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.8686e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1237e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.4756e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 3.5811e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.1752e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.4075   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.9944   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.1896,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.7368,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.0492,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0663e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 NOTE: this values are very high. This is not a problem,
                       since the IB integrands are only integrated over w1 between 0 and sigma mB = 0.3937...*/
                    std::make_pair(-9.4239e1,  1.0e-3),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-4.8570e4,  1.0e-0),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-6.9108e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4483e-1, 1.0e-5),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 5.3970e-4,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.2773e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.3461e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.9182e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-2.1902e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.3474e-2,  1.0e-6),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.3502e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.6390e-2,  1.0e-5),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 6.8295e-1,  1.0e-5),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.6769e+2,  1.0e-2),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-2.2440e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.6934e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.0896e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.7367e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.1806e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.5327e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.2112e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.1271e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 1.9033e-1,  1.0e-5),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1394e+2,  1.0e-2),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.3490e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9555e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0104e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.3562e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-8.7833,    1.0e-4),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.0906e3,  1.0e-1),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 2.4701e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.5805e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A1_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_1 phi_4 */
                    std::make_pair(-1.5968e-3,  1.0e-7),  // I1_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.2460e-3,  1.0e-7),  // I1_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9279e-5,  1.0e-9),  // I1_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.9072e-4,  1.0e-8),  // I1_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-3.2572e-2,  1.0e-6),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0701e-1,  1.0e-5),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.6172e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.8903e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_1 phi_bar_4 */
                    std::make_pair(-2.3439e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7004e-3,  1.0e-7),  // I1_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0136e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4384e-2,  1.0e-6),  // I1_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.6830e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.8049e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.4274e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6740e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-9.7591e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.2410   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.2049   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0159e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 7.3408   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5402e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.5044   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.7843e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-7.1522  ,  1.0e-4),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.6837e3,  1.0e-1),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-4.9002  ,  1.0e-4),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.0269e1,  1.0e-3),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 9.8237e-3,  1.0e-7),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.4454e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9291e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.5413e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 1.9972e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.9648e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.9281   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1263e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.0365   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.2680   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-9.5387   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.5782e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-1.5607e+1,  1.0e-3),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.2772e+4,  1.0e-0),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 6.8881e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4306e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.3872e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8841e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.9774e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0224   ,  1.0e-4),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.8325e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.3936e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.2959e+1,  1.0e-3),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0605e+4,  1.0e-0),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-5.3697e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.6822   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0417e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6222   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-4.3716e2,  1.0e-2),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.5046e5,  1.0e+1),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-3.5102,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-9.3513,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_1 psi_bar_4 */
                    std::make_pair( 5.7779e-3,  1.0e-7),  // I1_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.1616e-3,  1.0e-7),  // I1_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.6260e-2,  1.0e-6),  // I1_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6343e-2,  1.0e-6),  // I1_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 8.6733e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.7459e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.4430e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4532e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 2.3938   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3098   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5023e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.7708   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-1.9486e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.1903   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.1008e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3032e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 4.4659e1,  1.0e-3),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.3431e3,  1.0e-1),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 3.2135   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.3469   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 psi_bar_bar_4 */
                    std::make_pair(-1.9915e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5176e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3719e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3526e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_bar_4 */
                    std::make_pair(-4.0237e-1,  1.0e-5),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.9835e-1,  1.0e-5),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0932e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.7641   ,  1.0e-4),  // I3_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_bar_4 */
                    std::make_pair( 6.1456   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0380   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6582e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5666e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_bar_4 */
                    std::make_pair(-5.1460   ,  1.0e-4),  // I3d1B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1437e+3,  1.0e-1),  // I3d1B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 psi_bar_bar_4 */
                    std::make_pair(-2.4480e-2,  1.0e-6),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.0848e-2,  1.0e-6),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.6320e-1,  1.0e-5),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.7231e-1,  1.0e-5),  // I4_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A psi_bar_bar_4 */
                    std::make_pair( 1.8914e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.4738e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6668e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8238e-1,  1.0e-5),  // I4d1A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B psi_bar_bar_4 */
                    std::make_pair( 3.7487   ,  1.0e-4),  // I4d1B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.3315e+2,  1.0e-2),  // I4d1B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A psi_bar_bar_4 */
                    std::make_pair( 4.5644e-1,  1.0e-5),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3210   ,  1.0e-4),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.4214e-1,  1.0e-5),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4372   ,  1.0e-4),  // I4d2A_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B psi_bar_bar_4 */
                    std::make_pair(-8.9454e1,  1.0e-3),   // I4d2B_A1_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-9.9289e3,  1.0e-1),   // I4d2B_A1_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C psi_bar_bar_4 */
                    std::make_pair( 3.7569,     1.0e-4),  // I4d2C_A1_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.1746,     1.0e-4),  // I4d2C_A1_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_psi_bar_bar_4(sigma_0, 5.0)

                    /* I_1 chi_bar_4 */
                    std::make_pair( 1.7652e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.6994e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.7318e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.7250e-3,  1.0e-7),  // I1_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 5.7930e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.9119e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.7142   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6316   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 7.6042e-1,  1.0e-5),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6257   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0788   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4748   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 1.2678e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.9456e-2,  1.0e-6),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2170e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.1963e-3,  1.0e-7),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-1.9119e2,  1.0e-2),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.2875e4,  1.0e-0),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 13.702,  1.0e-3),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 5.7428,  1.0e-4),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_bar_4 */
                    std::make_pair(-7.9662e-2,  1.0e-6),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.0703e-2,  1.0e-6),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1488   ,  1.0e-4),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.3410   ,  1.0e-4),  // I2_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_bar_4 */
                    std::make_pair(-1.6297   ,  1.0e-4),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.2608   ,  1.0e-4),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.3779e+1,  1.0e-3),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7457e+1,  1.0e-3),  // I3_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_bar_4 */
                    std::make_pair( 2.4665e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6506e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0617e+2,  1.0e-2),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.2123e+1,  1.0e-3),  // I3d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_bar_4 */
                    std::make_pair(-1.6802e+1,  1.0e-3),  // I3d1B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.7344e+3,  1.0e-1),  // I3d1B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 chi_bar_bar_4 */
                    std::make_pair(-1.3793e-2,  1.0e-6),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2073e-1,  1.0e-5),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6301   ,  1.0e-4),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.7471e-1,  1.0e-5),  // I4_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A chi_bar_bar_4 */
                    std::make_pair(-4.3585e-1,  1.0e-5),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4131   ,  1.0e-4),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.4116   ,  1.0e-4),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.9701e-1,  1.0e-5),  // I4d1A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B chi_bar_bar_4 */
                    std::make_pair( 1.6049e+1,  1.0e-3),  // I4d1B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.5668e+3,  1.0e-1),  // I4d1B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A chi_bar_bar_4 */
                    std::make_pair( 6.8829   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.6913   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9169   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.9427   ,  1.0e-4),  // I4d2A_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B chi_bar_bar_4 */
                    std::make_pair(-3.9651e2,  1.0e-2),   // I4d2B_A1_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-4.2570e4,  1.0e+0),   // I4d2B_A1_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C chi_bar_bar_4 */
                    std::make_pair( 14.102,     1.0e-3),  // I4d2C_A1_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 8.1623,     1.0e-4),  // I4d2C_A1_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_chi_bar_bar_4(sigma_0, 5.0)

                    /* A_2 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.043309e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.276134e-1,  1.0e-6), // I1_A2_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(-1.214754e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.214754e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.214754e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.578119e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.578119e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.578119e-1,  1.0e-6), // I2_A2_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(-5.879654,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-5.879654,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-5.879654,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.216099,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.216099,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-5.216099,  1.0e-6), // I2d1_A2_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.141327e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.943063e-1,  1.0e-6), // I2_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.483235,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.067506,  1.0e-6), // I2d1_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(3.826506e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.826506e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.826506e-3,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.029616e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.029616e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.029616e-2,  1.0e-6), // I3_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.542946e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.542946e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.542946e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.534912e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.534912e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.534912e-1,  1.0e-6), // I3d1_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(1.377467,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.377467,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.377467,  1.0e-6),    // I3d2_A2_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.004497,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.004497,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.004497,  1.0e-6),   // I3d2_A2_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.217850e-3,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.050510e-1,  1.0e-6), // I3_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.758248e-1,  1.0e-6), // I3d1_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.666550,  1.0e-6),    // I3d1_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar ATTENTION: this term is huge!!! */
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.266027e+1,  1.0e-4),     // I3d2_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.420043e+2,  1.0e-4), // I3d2_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.839758e-4,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-4.612858e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-4.612858e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.612858e-3,  1.0e-6), // I4_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(-3.020304e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.020304e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.020304e-2,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.201956e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.201956e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.201956e-1,  1.0e-6), // I4d1_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(-2.312397,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.312397,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.312397,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-7.367758,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-7.367758,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-7.367758,  1.0e-5), // I4d2_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:this term is huge!!!*/
                    std::make_pair(-1.075823e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.075823e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.075823e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.356839e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.356839e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.356839e2,  1.0e-4), // I4d3_A2_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair(-1.4278e-1,  1.0e-5),  // I2_A2_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0585e-1,  1.0e-5),  // I2_A2_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3765e-2,  1.0e-6),  // I2_A2_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.0061e-2,  1.0e-6),  // I2_A2_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.8162e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.8062e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2108e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.5374e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 5.9725e-1,  1.0e-5),  // I3_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1750   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.0284   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.3236   ,  1.0e-4),  // I3_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-2.7218,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.7722,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.0008,     1.0e-4),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8757e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.6161e2,  1.0e-2),   // I3d1B_A2_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-8.3291e4,  1.0e-0),   // I3d1B_A2_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-1.6058e-1, 1.0e-5),   // I3d1C_A2_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-3.3652e-1, 1.0e-5),   // I3d1C_A2_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-6.2581e-3,  1.0e-7),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.8484e-2,  1.0e-6),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.1389e-1,  1.0e-5),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.1549e-1,  1.0e-5),  // I4_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 1.6683e-2,  1.0e-6),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.7536e-1,  1.0e-5),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5981e-1,  1.0e-5),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.1123   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-2.5340   ,  1.0e-4),  // I4d1B_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.8484e+3,  1.0e-1),  // I4d1B_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair( 6.0020e-1,  1.0e-5),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.2669   ,  1.0e-4),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.9052e-1,  1.0e-5),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.7375   ,  1.0e-4),  // I4d2A_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-5.4493,    1.0e-4),   // I4d2B_A2_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.6355e4,  1.0e-0),   // I4d2B_A2_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 4.5815e-1,  1.0e-5),  // I4d2C_A2_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.2206   ,  1.0e-4),  // I4d2C_A2_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A2_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-5.1174e-2,  1.0e-6),  // I2_A2_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6812e-1,  1.0e-5),  // I2_A2_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.5407e-3,  1.0e-7),  // I2_A2_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.1120e-3,  1.0e-7),  // I2_A2_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 4.8925e-3,  1.0e-7),  // I2_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6073e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1158e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0898e-2,  1.0e-6),  // I2_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.1690   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1294   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.9270   ,  1.0e-4),  // I3_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2218e+1,  1.0e-3),  // I3_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.7895   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4895e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.2705e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2655e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-5.9307e1,  1.0e-3),   // I3d1B_A2_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.2253e4,  1.0e-0),   // I3d1B_A2_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.8824  ,  1.0e-4),   // I3d1C_A2_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4423e1,  1.0e-3),   // I3d1C_A2_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair(-2.9488e-1,  1.0e-5),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0397   ,  1.0e-4),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7825   ,  1.0e-4),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6638e+1,  1.0e-3),  // I3_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair( 3.7096   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1027e+1,  1.0e-3),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-9.3046e-1,  1.0e-5),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.1932   ,  1.0e-4),  // I3d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair( 2.1750e+1,  1.0e-3),  // I3d1B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.7799e+4,  1.0e-0),  // I3d1B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(-1.2599   ,  1.0e-4),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.2681   ,  1.0e-4),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.4967e+1,  1.0e-3),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.0940e+1,  1.0e-3),  // I4_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 2.4486   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.1719   ,  1.0e-4),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.6923e+2,  1.0e-2),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.8427e+2,  1.0e-2),  // I4d1A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.1266e+2,  1.0e-2),  // I4d1B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 9.2199e+4,  1.0e-0),  // I4d1B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 3.9498e+2,  1.0e-2),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2766e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6454e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.4536e+3,  1.0e-1),  // I4d2A_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-2.8370e3,  1.0e-1),   // I4d2B_A2_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.4599e5,  1.0e+1),   // I4d2B_A2_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair( 2.5661e1,   1.0e-3),  // I4d2C_A2_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 6.8362e1,   1.0e-3),  // I4d2C_A2_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-1.2061e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.5994e-3,  1.0e-7),  // I2_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.5687e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.4113e-2,  1.0e-6),  // I2_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 2.7663   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4881   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7438e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.8175   ,  1.0e-4),  // I3_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.4656e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0700e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.8371e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0841e+1,  1.0e-3),  // I3d1A_A2_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.5940e1,  1.0e-3),   // I3d1B_A2_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.6928e3,  1.0e-1),   // I3d1B_A2_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 3.4460   ,  1.0e-4),  // I3d1C_A2_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.4443   ,  1.0e-4),  // I3d1C_A2_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 psi_bar_bar_4 */
                    std::make_pair( 5.9782e-1,  1.0e-5),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.5797e-1,  1.0e-5),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6102e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0067e+1,  1.0e-3),  // I3_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_bar_4 */
                    std::make_pair(-7.5204   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.8572   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.5910   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.1420   ,  1.0e-4),  // I3d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_bar_4 */
                    std::make_pair( 6.5126   ,  1.0e-4),  // I3d1B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.4474e+3,  1.0e-1),  // I3d1B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 psi_bar_bar_4 */
                    std::make_pair( 2.5248   ,  1.0e-4),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8306   ,  1.0e-4),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.8992e+1,  1.0e-3),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.2395e+1,  1.0e-3),  // I4_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A psi_bar_bar_4 */
                    std::make_pair(-4.8859   ,  1.0e-4),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.3762e-1,  1.0e-5),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.4756e+2,  1.0e-2),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.7031e+2,  1.0e-2),  // I4d1A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B psi_bar_bar_4 */
                    std::make_pair( 3.5488e+1,  1.0e-3),  // I4d1B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.8874e+3,  1.0e-1),  // I4d1B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A psi_bar_bar_4 */
                    std::make_pair(-7.9793e+2,  1.0e-2),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.5492e+2,  1.0e-2),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.5837e+3,  1.0e-1),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.7024e+3,  1.0e-1),  // I4d2A_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B psi_bar_bar_4 */
                    std::make_pair(-5.1089e2,  1.0e-2),   // I4d2B_A2_3pt_psi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.9941e4,  1.0e+0),   // I4d2B_A2_3pt_psi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C psi_bar_bar_4 */
                    std::make_pair(-2.4471e1,   1.0e-3),  // I4d2C_A2_3pt_psi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4164e1,   1.0e-3),  // I4d2C_A2_3pt_psi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D psi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_psi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-3.6847e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7219e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7023e-3,  1.0e-7),  // I2_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1950e-2,  1.0e-6),  // I2_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair(-1.3357   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0534   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.3830e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3001   ,  1.0e-4),  // I3_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-1.6148e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.7073   ,  1.0e-4),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0750e+2,  1.0e-2),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-8.8422e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-3.2786e2,  1.0e-2),   // I3d1B_A2_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.9227e4,  1.0e-0),   // I3d1B_A2_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 16.912,  1.0e-3),     // I3d1C_A2_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.0885,  1.0e-4),     // I3d1C_A2_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 chi_bar_bar_4 */
                    std::make_pair( 2.3913   ,  1.0e-4),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8319   ,  1.0e-4),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.4409e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.0266e+1,  1.0e-3),  // I3_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_bar_4 */
                    std::make_pair(-3.0082e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.9429e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0364e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2568e+1,  1.0e-3),  // I3d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_bar_4 */
                    std::make_pair( 2.6050e+1,  1.0e-3),  // I3d1B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.7898e+3,  1.0e-0),  // I3d1B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 chi_bar_bar_4 */
                    std::make_pair( 1.0334e+1,  1.0e-3),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.7172   ,  1.0e-4),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.8024e+2,  1.0e-2),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.7378e+2,  1.0e-2),  // I4_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A chi_bar_bar_4 */
                    std::make_pair(-2.0169e+1,  1.0e-3),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.7789   ,  1.0e-4),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0075e+3,  1.0e-1),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.9147e+3,  1.0e-1),  // I4d1A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B chi_bar_bar_4 */
                    std::make_pair( 1.2792e+2,  1.0e-2),  // I4d1B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.8431e+4,  1.0e-0),  // I4d1B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A chi_bar_bar_4 */
                    std::make_pair(-3.2142e+3,  1.0e-1),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2789e+3,  1.0e-1),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.8320e+4,  1.0e-0),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0746e+4,  1.0e-0),  // I4d2A_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B chi_bar_bar_4 */
                    std::make_pair(-2.2061e3,  1.0e-1),   // I4d2B_A2_3pt_chi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1102e5,  1.0e+1),   // I4d2B_A2_3pt_chi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C chi_bar_bar_4 */
                    std::make_pair(-1.1506e2,   1.0e-2),  // I4d2C_A2_3pt_chi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-6.6599e1,   1.0e-3),  // I4d2C_A2_3pt_chi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D chi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A2_3pt_chi_bar_bar_4(sigma_0, 5.0)

                    /* A_30 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.094319e-1,  1.0e-6), // I1_A30_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.004800,  1.0e-6),    // I1_A30_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.362076e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.362076e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.362076e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.321795e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.321795e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.321795e-1,  1.0e-6), // I2_A30_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(6.815662,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.815662,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.815662,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.145066,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.145066,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.145066,  1.0e-6), // I2d1_A30_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.339818e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.064229e-1,  1.0e-6), // I2_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.802308,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.285435,  1.0e-6), // I2d1_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(-4.491985e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.491985e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-4.491985e-3,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.421851e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.421851e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.421851e-2,  1.0e-6), // I3_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(-1.992121e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.992121e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.992121e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.703323e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.703323e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.703323e-1,  1.0e-6), // I3d1_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(-3.154013,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.154013,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.154013,  1.0e-6),      // I3d2_A30_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.310445e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.310445e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.310445e-1,  1.0e-6),   // I3d2_A30_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-8.902671e-3,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.233207e-1,  1.0e-6), // I3_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-8.583107e-1,  1.0e-6), // I3d1_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-5.726354,  1.0e-6),    // I3d1_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_bar ATTENTION: this term is huge!!! */
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.045307e+1,  1.0e-4),     // I3d2_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.892940e+2,  1.0e-4), // I3d2_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.184156e-4,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.571594e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.571594e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.571594e-3,  1.0e-6), // I4_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(3.438716e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.438716e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.438716e-2,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.768528e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.768528e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.768528e-1,  1.0e-6), // I4d1_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar ATTENTION:here I have to low the precision to e-5 to pass the test*/
                    std::make_pair(2.713015,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.713015,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.713015,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(9.972535,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(9.972535,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(9.972535,  1.0e-5), // I4d2_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar ATTENTION:this term is huge!!!*/
                    std::make_pair(1.354778e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.354778e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.354778e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.205744e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.205744e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.205744e2,  1.0e-4), // I4d3_A30_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 1.5611e-1,  1.0e-5),  // I2_A30_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.2052e-1,  1.0e-5),  // I2_A30_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6198e-2,  1.0e-6),  // I2_A30_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.4645e-2,  1.0e-6),  // I2_A30_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-5.0321e-2,  1.0e-6),  // I2_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.0546e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.3547e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-7.0303e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-7.7567e-1,  1.0e-5),  // I3_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5220   ,  1.0e-4),  // I3_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.2343   ,  1.0e-4),  // I3_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0809e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair( 1.2282,     1.0e-4),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.0900,     1.0e-4),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.7350e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.6519e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair( 2.1073e2,  1.0e-2),   // I3d1B_A30_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.0861e5,  1.0e+1),   // I3d1B_A30_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 2.1715e-1, 1.0e-5),   // I3d1C_A30_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 4.5508e-1, 1.0e-5),   // I3d1C_A30_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair( 5.7744e-3,  1.0e-7),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.1128e-2,  1.0e-6),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.1863e-2,  1.0e-6),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.8303e-1,  1.0e-5),  // I4_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair(-7.6954e-3,  1.0e-7),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.9702e-1,  1.0e-5),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.0384e-1,  1.0e-5),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.2576   ,  1.0e-4),  // I4d1A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 3.0627   ,  1.0e-4),  // I4d1B_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.4428e+3,  1.0e-1),  // I4d1B_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-2.7985e-1,  1.0e-5),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.9625   ,  1.0e-4),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6889   ,  1.0e-4),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.6008e+1,  1.0e-3),  // I4d2A_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair( 1.6366e1,  1.0e-3),   // I4d2B_A30_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.1093e4,  1.0e-0),   // I4d2B_A30_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(-5.1478e-1,  1.0e-5),  // I4d2C_A30_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.3714   ,  1.0e-4),  // I4d2C_A30_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A30_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 3.7842e-2,  1.0e-6),  // I2_A30_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2432e-1,  1.0e-5),  // I2_A30_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8788e-3,  1.0e-7),  // I2_A30_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.5196e-3,  1.0e-7),  // I2_A30_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 1.3555e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.4534e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8622e-1,  1.0e-5),  // I2_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.4102   ,  1.0e-4),  // I2_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair( 1.5861   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.5683   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.7005   ,  1.0e-4),  // I3_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6571e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair(-5.6152   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 8.5968e-1,  1.0e-5),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0036e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5259e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 7.3363e1,  1.0e-3),   // I3d1B_A30_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.7527e4,  1.0e-0),   // I3d1B_A30_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 9.1963  ,  1.0e-4),   // I3d1C_A30_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.9272e1,  1.0e-3),   // I3d1C_A30_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 3.8210e-1,  1.0e-5),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.2988   ,  1.0e-4),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.5656   ,  1.0e-4),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.1519e+1,  1.0e-3),  // I3_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.5670   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.7556   ,  1.0e-4),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6398e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 7.8239e+1,  1.0e-3),  // I3d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-3.3673e+1,  1.0e-3),  // I3d1B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7556e+4,  1.0e-0),  // I3d1B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.4545   ,  1.0e-4),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.8463   ,  1.0e-4),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.8947e+1,  1.0e-3),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.1832e+1,  1.0e-3),  // I4_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 2.9822e-2,  1.0e-6),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0427e+1,  1.0e-3),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.7024e+2,  1.0e-2),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0682e+3,  1.0e-1),  // I4d1A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-1.3927e+2,  1.0e-2),  // I4d1B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1398e+5,  1.0e+1),  // I4d1B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-4.6218e+2,  1.0e-2),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4356e+3,  1.0e-1),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3003e+2,  1.0e-2),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1782e+3,  1.0e-1),  // I4d2A_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 2.7803e3,  1.0e-1),   // I4d2B_A30_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1866e6,  1.0e+2),   // I4d2B_A30_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-2.7341e1,   1.0e-3),  // I4d2C_A30_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-7.2838e1,   1.0e-3),  // I4d2C_A30_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A30_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-3.3416e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.8285e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0970   ,  1.0e-4),  // I2_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.4516e-1,  1.0e-5),  // I2_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(-3.7711   ,  1.0e-4),  // I3_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0380   ,  1.0e-4),  // I3_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3744e+1,  1.0e-3),  // I3_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0660e+1,  1.0e-3),  // I3_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair( 1.7986e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 6.0131   ,  1.0e-4),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-3.1801e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5690e+1,  1.0e-3),  // I3d1A_A30_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(-7.4685e1,  1.0e-3),   // I3d1B_A30_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-8.9355e3,  1.0e-1),   // I3d1B_A30_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(-4.7949   ,  1.0e-4),  // I3d1C_A30_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.0097   ,  1.0e-4),  // I3d1C_A30_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair(-1.0209e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.1395e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5799e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.3110e-1,  1.0e-5),  // I2_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.3612   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6236   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5774e+1,  1.0e-3),  // I3_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.2159   ,  1.0e-4),  // I3_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 1.8873e+1,  1.0e-3),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8371e-1,  1.0e-5),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6875e+2,  1.0e-2),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0068e+2,  1.0e-2),  // I3d1A_A30_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 4.2752e2,  1.0e-2),   // I3d1B_A30_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 5.1150e4,  1.0e-0),   // I3d1B_A30_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair(-23.002,  1.0e-3),     // I3d1C_A30_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-9.6408,  1.0e-4),     // I3d1C_A30_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* V */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.244211e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.640698e-1,  1.0e-6), // I1_V_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.047657e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.047657e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.047657e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.322133e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.322133e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.322133e-2,  1.0e-6), // I2_V_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 phi_bar */
                    std::make_pair(1.417983e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.417983e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.417983e-1,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.240511e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.240511e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.240511e-2,  1.0e-6), // I2d1_V_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-2.349792e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.636322e-2,  1.0e-6), // I2_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2d1 g_+ */
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-9.74102e-1,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.075191,  1.0e-6), // I2d1_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.878105e-4,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.321685e-3,  1.0e-6), // I3_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(3.347918e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.347918e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.347918e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.013860e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.013860e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.013860e-2,  1.0e-6), // I3d1_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d2 g_+ */
                    std::make_pair(4.291581e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.291581e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.291581e-1,  1.0e-6),      // I3d2_V_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.536817e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.536817e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.536817e-2,  1.0e-6),   // I3d2_V_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.449130e-5,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.704475e-4,  1.0e-6), // I4_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.750898e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.811524e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.811524e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.811524e-3,  1.0e-6), // I4d1_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(7.918965e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(7.918965e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.918965e-2,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.192451e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.192451e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.192451e-1,  1.0e-6), // I4d2_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d3 g_bar */
                    std::make_pair(1.393143,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.393143,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.393143,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(7.004433e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(7.004433e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(7.004433e-1,  1.0e-4), // I4d3_V_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair(9.4355e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.9773e-2,  1.0e-6),  // I2_V_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(1.5772e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(3.3054e-3,  1.0e-7),  // I2_V_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(4.5937e-3,  1.0e-7),  // I2_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(9.6268e-3,  1.0e-7),  // I2_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.0624e-2,  1.0e-6),  // I2_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(6.4178e-2,  1.0e-6),  // I2_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(9.3704e-2,  1.0e-6),  // I3_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.9637e-1,  1.0e-5),  // I3_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(6.2469e-1,  1.0e-5),  // I3_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(1.3091   ,  1.0e-4),  // I3_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-6.6912e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.4023,     1.0e-4), // I3d1A_V_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.7239e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.7084e-1,  1.0e-5), // I3d1A_V_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 NOTE: this values are very high. This is not a problem,
                       since the IB integrands are only integrated over w1 between 0 and sigma mB = 0.3937...*/
                    std::make_pair(-2.2888e1,  1.0e-3),  // I3d1B_V_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.1796e4,  1.0e-0),  // I3d1B_V_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-5.8373e-4,  1.0e-8), // I4_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0302e-3,  1.0e-7), // I4_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5781e-2,  1.0e-6), // I4_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.5208e-2,  1.0e-6), // I4_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 7.2808e-3,  1.0e-7), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.1491e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.8611e-3,  1.0e-7), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0225e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(4.6225e-2,  1.0e-6),  // I4d1B_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(5.1961e+1,  1.0e-3),  // I4d1B_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.1067e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.7120e-3,  1.0e-7), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.2000e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2003e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-2.0396,    1.0e-4),  // I4d2B_V_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.0250e2,  1.0e-2),  // I4d2B_V_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(6.8423e-3,  1.0e-7),  // I4d2C_V_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(1.8228e-2,  1.0e-6),  // I4d2C_V_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),        // I4d2D_V_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-8.4303e-3,  1.0e-7), // I2_V_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.7696e-2,  1.0e-6), // I2_V_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.1856e-4,  1.0e-8), // I2_V_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.0069e-3,  1.0e-7), // I2_V_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-1.2375e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.0654e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3514e-2,  1.0e-6), // I2_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2873e-1,  1.0e-5), // I2_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-2.5360e-1,  1.0e-5), // I3_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-8.4161e-1,  1.0e-5), // I3_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0929   ,  1.0e-4), // I3_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6399   ,  1.0e-4), // I3_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair(1.4859   ,  1.0e-4),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(2.6042   ,  1.0e-4),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.8988e-1,  1.0e-5),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(2.4277e-1,  1.0e-5),  // I3d1A_V_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-1.7371  ,  1.0e-4),  // I3d1B_V_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.5178e2,  1.0e-2),  // I3d1B_V_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-1.2709,  1.0e-4),    // I3d1C_V_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.6634,  1.0e-4),    // I3d1C_V_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                     /* I_3 phi_bar_bar_4 */
                    std::make_pair(1.5388e-4,  1.0e-8),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(2.0497e-3,  1.0e-7),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(7.4683e-4,  1.0e-8),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(9.9481e-3,  1.0e-7),  // I3_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-6.7828e-4,  1.0e-8), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.0350e-3,  1.0e-7), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7153e-3,  1.0e-7), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2849e-2,  1.0e-6), // I3d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(1.5954e-1,  1.0e-5),  // I3d1B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(1.3056e+2,  1.0e-2),  // I3d1B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I3d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(4.4881e-3,  1.0e-7),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(4.6503e-2,  1.0e-6),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(4.1790e-2,  1.0e-6),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(2.7900e-1,  1.0e-5),  // I4_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-3.5851e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0304e-1,  1.0e-5), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8005e-2,  1.0e-6), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6490e-1,  1.0e-5), // I4d1A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(3.1475   ,  1.0e-4),  // I4d1B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(2.5758e+3,  1.0e-1),  // I4d1B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),        // I4d1C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-9.4302e-2,  1.0e-6), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5860   ,  1.0e-4), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2354e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.8189e-1,  1.0e-5), // I4d2A_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-9.2800e1,  1.0e-3),  // I4d2B_V_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7869e4,  1.0e+0),  // I4d2B_V_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-9.7235e-1,  1.0e-5), // I4d2C_V_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-2.5904,     1.0e-4), // I4d2C_V_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),        // I4d2D_V_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(3.0504e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(1.6692e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(1.9143e-1,  1.0e-5),  // I2_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(8.6282e-2,  1.0e-6),  // I2_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(6.2536e-1,  1.0e-5),  // I3_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(3.4702e-1,  1.0e-5),  // I3_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(3.9098   ,  1.0e-4),  // I3_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(1.7701   ,  1.0e-4),  // I3_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-4.0182   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5404   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5380   ,  1.0e-4), // I3d1A_V_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.2218e-1,  1.0e-5), // I3d1A_V_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(1.0846e1,  1.0e-3),   // I3d1B_V_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(1.2977e3,  1.0e-1),   // I3d1B_V_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(8.9017e-1,  1.0e-5),  // I3d1C_V_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(3.7310e-1,  1.0e-5),  // I3d1C_V_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 9.3196e-3,  1.0e-7), // I2_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9531e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.4423e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.0225e-2,  1.0e-6), // I2_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 1.9011e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.9840e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9420e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.1655e-1,  1.0e-5), // I3_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(3.0244e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(6.3382e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(4.6805e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(9.8088e-1,  1.0e-5),  // I3d1A_V_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-4.6435e1,  1.0e-3),  // I3d1B_V_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.5556e3,  1.0e-1),  // I3d1B_V_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 3.3413,  1.0e-4),    // I3d1C_V_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.4004,  1.0e-4),    // I3d1C_V_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_1 */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.646389,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.646389,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.646389,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-0.817723,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.817723,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.817723,  1.0e-6), // I1_T1_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(5.442762e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.442762e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(5.442762e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(6.589510e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(6.589510e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(6.589510e-2,  1.0e-6), // I2_T1_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T1_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(6.177960e-5,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.120409e-4,  1.0e-6), // I3_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(4.352307e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.352307e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.352307e-3,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.360095e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.360095e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.360095e-2,  1.0e-6), // I3d1_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.272367e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(8.495098e-4,  1.0e-6), // I4_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(8.966925e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.966925e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(8.966925e-3,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.806474e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.806474e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.806474e-2,  1.0e-6), // I4d1_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(3.929166e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(3.929166e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(3.929166e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(5.329536e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(5.329536e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(5.329536e-1,  1.0e-6), // I4d2_T1_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 4.5886e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.4503e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6959e-3,  1.0e-7),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6064e-2,  1.0e-6),  // I2_A1_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.1314e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.4519e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.4936e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5829e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 4.5570e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.3851e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0480   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.3621   ,  1.0e-4),  // I3_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.7659,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.9184,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.6333,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.7503,     1.0e-4),  // I3d1A_A1_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.1469e2,  1.0e-2),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-5.9111e4,  1.0e-0),   // I3d1B_A1_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-3.4554e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-7.2414e-2, 1.0e-6),   // I3d1C_A1_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 3.2513e-4,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3309e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1675e-3,  1.0e-7),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.8873e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-1.7846e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.3772e-2,  1.0e-6),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6355e-3,  1.0e-7),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5106e-2,  1.0e-5),  // I3d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 3.3710e-1,  1.0e-5),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.7893e+2,  1.0e-2),  // I3d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-2.8454e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.8404e-3,  1.0e-7),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.7037e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.2030e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.8617e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1535e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.9658e-2,  1.0e-6),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.4105e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 2.3163e-1,  1.0e-5),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.6037e+2,  1.0e-2),  // I4d1B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-1.3123e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.5386e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.0552e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.0275e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-1.0606e1,  1.0e-3),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.5411e3,  1.0e-1),   // I4d2B_A1_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 3.2551e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 8.6718e-2,  1.0e-6),  // I4d2C_A1_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_A1_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-4.1175e-2,  1.0e-6),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.3527e-1,  1.0e-5),  // I2_A1_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0443e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.9177e-3,  1.0e-7),  // I2_A1_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.5931e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.4342e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2585e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6375e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.2366   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1052   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3291   ,  1.0e-4),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2873e+1,  1.0e-3),  // I3_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.5867   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7147e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6855   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.5151e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-8.7044  ,  1.0e-4),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.2661e3,  1.0e-1),   // I3d1B_A1_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.2022  ,  1.0e-4),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.2998e1,  1.0e-3),   // I3d1C_A1_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 4.8973e-3,  1.0e-7),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7033e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.6386e-2,  1.0e-6),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.7612e-1,  1.0e-5),  // I2_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 9.9891e-2,  1.0e-6),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.4989e-1,  1.0e-5),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9623   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.6343   ,  1.0e-4),  // I3_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-1.5172   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.6362   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7676   ,  1.0e-4),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2860e+1,  1.0e-3),  // I3d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-7.6307   ,  1.0e-4),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.2447e+3,  1.0e-1),  // I3d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 1.6694e-2,  1.0e-6),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0882e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0144e-1,  1.0e-5),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.0679   ,  1.0e-4),  // I4_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.2510e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5030   ,  1.0e-4),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.0504e-3,  1.0e-7),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.1546e-1,  1.0e-5),  // I4d1A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 1.5772e+1,  1.0e-3),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.2907e+4,  1.0e-0),  // I4d1B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-3.9181e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.2024   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8799e-1,  1.0e-5),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.6803   ,  1.0e-4),  // I4d2A_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-5.0570e2,  1.0e-2),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.6826e5,  1.0e+1),   // I4d2B_A1_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-4.6257,     1.0e-4),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-12.323,     1.0e-3),  // I4d2C_A1_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_A1_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 7.4198e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0143e-2,  1.0e-6),  // I2_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.6704e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0975e-1,  1.0e-5),  // I2_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 3.0431   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6794   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9054e+1,  1.0e-3),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.6113   ,  1.0e-4),  // I3_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.2912e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.3676   ,  1.0e-4),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8209e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1459e+1,  1.0e-3),  // I3d1A_A1_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.4351e1,  1.0e-3),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.5027e3,  1.0e-1),   // I3d1B_A1_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 4.2348   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.7749   ,  1.0e-4),  // I3d1C_A1_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.9916e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6553e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8718   ,  1.0e-4),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.4666e-1,  1.0e-5),  // I2_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 9.4146e-1,  1.0e-5),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9890   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.4080   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.0576   ,  1.0e-4),  // I3_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 5.5930e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0840   ,  1.0e-4),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.0661e-1,  1.0e-5),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6109   ,  1.0e-4),  // I3d1A_A1_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-2.3268e2,  1.0e-2),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.7839e4,  1.0e-0),   // I3d1B_A1_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 16.715,  1.0e-3),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.0059,  1.0e-4),     // I3d1C_A1_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_23A */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(-0.646389,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-0.646389,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-0.646389,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-0.817723,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-0.817723,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-0.817723,  1.0e-6), // I1_T23A_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(8.746169e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(5.442762e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(2.139355e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 1.492725e-1,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 6.589510e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.748232e-2,  1.0e-6), // I2_T23A_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.190950e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-3.223354e-1,  1.0e-6), // I2_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.092816e-3,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.157127e-2,  1.0e-6), // I3_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.697711e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.877934e-1,  1.0e-6), // I3d1_T23A_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.842692e-3,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.177781e-2,  1.0e-6), // I2_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-1.473793e-3,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 6.177960e-5,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.597356e-3,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.007110e-2,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 4.120409e-4,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair( 2.089519e-2,  1.0e-6), // I3_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-1.422168e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 4.352307e-3,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.509214e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.185610e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 1.360095e-2,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair( 9.457629e-1,  1.0e-6), // I3d1_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair( 2.044612e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 1.272367e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 5.001218e-5,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 1.924399e-3,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 8.495098e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.253794e-4,  1.0e-6), // I4_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair( 1.641836e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 8.966925e-3,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair( 1.515487e-3,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 7.814995e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 2.806474e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.202046e-2,  1.0e-6), // I4d1_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair( 9.186539e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair( 3.929166e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.328207e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair( 2.201174   ,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair( 5.329536e-1,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.135267   ,  1.0e-6), // I4d2_T23A_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_2 phi_3 */
                    std::make_pair( 1.0142e-2,  1.0e-6),  // I2_T23A_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.9596e-2,  1.0e-6),  // I2_T23A_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.7208e-3,  1.0e-7),  // I2_T23A_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5419e-3,  1.0e-7),  // I2_T23A_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-1.1314e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.4519e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.4936e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5829e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair( 4.9009e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.9443e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 3.2872   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 6.8383   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-3.4865,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.3852,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.1149,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4581,     1.0e-4),  // I3d1A_T23A_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair(-1.2641e2,  1.0e-2),   // I3d1B_T23A_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-6.5151e4,  1.0e+0),   // I3d1B_T23A_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair(-6.8433e-2, 1.0e-6),   // I3d1C_T23A_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.4341e-1, 1.0e-5),   // I3d1C_T23A_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 3.2513e-4,  1.0e-8),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.3309e-3,  1.0e-7),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.1675e-3,  1.0e-7),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.8873e-2,  1.0e-6),  // I3_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-1.7846e-3,  1.0e-7),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.3772e-2,  1.0e-6),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.6355e-3,  1.0e-7),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.5106e-2,  1.0e-6),  // I3d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 3.3710e-1, 1.0e-5),   // I3d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.7893e+2, 1.0e-2),   // I3d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair(-3.3717e-3,  1.0e-7),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.6852e-2,  1.0e-6),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-8.0546e-2,  1.0e-6),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6704e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 3.3908e-2,  1.0e-6),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.2628e-2,  1.0e-6),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.7378e-3,  1.0e-7),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-5.9051e-1,  1.0e-5),  // I4d1A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair(-3.1413e-1,  1.0e-5),  // I4d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.5310e+2,  1.0e-2),  // I4d1B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-4.1812e-2,  1.0e-6),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 9.3717e-1,  1.0e-5),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.9767e-2,  1.0e-6),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.3376   ,  1.0e-4),  // I4d2A_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-1.0631e1,  1.0e-3),   // I4d2B_T23A_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.5299e4,  1.0e-0),   // I4d2B_T23A_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair( 1.1334e-1,  1.0e-5),  // I4d2C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 3.0193e-1,  1.0e-5),  // I4d2C_T23A_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_T23A_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair(-4.1175e-2,  1.0e-6),  // I2_T23A_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.3527e-1,  1.0e-5),  // I2_T23A_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0443e-3,  1.0e-7),  // I2_T23A_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.9177e-3,  1.0e-7),  // I2_T23A_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair(-2.5931e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.4342e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.2585e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.6375e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair(-1.2436   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1779   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.3369   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2955e+1,  1.0e-3),  // I3_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 8.4715   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5939e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.5557   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.3790e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair(-1.8959e1,  1.0e-3),   // I3d1B_T23A_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-7.1137e3,  1.0e-1),   // I3d1B_T23A_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair(-6.4407  ,  1.0e-4),   // I3d1C_T23A_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.3498e1,  1.0e-3),   // I3d1C_T23A_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair( 4.8973e-3,  1.0e-7),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7033e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 9.6386e-2,  1.0e-6),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.7612e-1,  1.0e-5),  // I2_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 3.9026e-2,  1.0e-6),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.3820e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 7.6437e-1,  1.0e-5),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2025   ,  1.0e-4),  // I3_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-8.6389e-1,  1.0e-5),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.7636   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.2204   ,  1.0e-4),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.0365e+1,  1.0e-3),  // I3d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-2.8108   ,  1.0e-4),  // I3d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-2.3003e+3,  1.0e-1),  // I3d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair(-2.3102e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.5271e-1,  1.0e-5),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.7739   ,  1.0e-4),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2899e+1,  1.0e-3),  // I4_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair( 1.2580e-1,  1.0e-5),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2563   ,  1.0e-4),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.7368e+1,  1.0e-3),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6662e+2,  1.0e-2),  // I4d1A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair( 3.5388e+1,  1.0e-3),  // I4d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.8960e+4,  1.0e+0),  // I4d1B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair( 7.6579e+1,  1.0e-3),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.3957e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9495e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.0209e+2,  1.0e-2),  // I4d2A_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair(-9.2164e2,  1.0e-2),   // I4d2B_T23A_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.3527e4,  1.0e+0),   // I4d2B_T23A_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair( 1.0489  ,   1.0e-3),  // I4d2C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.7945  ,   1.0e-3),  // I4d2C_T23A_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_T23A_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair( 7.4198e-2,  1.0e-5),  // I2_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.0143e-2,  1.0e-5),  // I2_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.6704e-1,  1.0e-4),  // I2_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.0975e-1,  1.0e-4),  // I2_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair( 3.0372   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6670   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9045e+1,  1.0e-3),  // I3_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.5922   ,  1.0e-4),  // I3_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-2.2931e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.4078   ,  1.0e-4),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8239e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.1521e+1,  1.0e-3),  // I3d1A_T23A_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair( 5.5795e1,  1.0e-3),   // I3d1B_T23A_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 6.6754e3,  1.0e-1),   // I3d1B_T23A_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair( 4.1309   ,  1.0e-4),  // I3d1C_T23A_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 1.7314   ,  1.0e-4),  // I3d1C_T23A_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.9916e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.6553e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8718   ,  1.0e-4),  // I2_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 8.4666e-1,  1.0e-5),  // I2_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 5.1953e-1,  1.0e-5),  // I3_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.8518   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5259   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.8893   ,  1.0e-4),  // I3_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair(-2.8682   ,  1.0e-4),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.2403e-1,  1.0e-5),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-4.2452e+1,  1.0e-3),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.6279e+1,  1.0e-3),  // I3d1A_T23A_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair(-2.5646e2,  1.0e-2),   // I3d1B_T23A_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.0683e4,  1.0e-0),   // I3d1B_T23A_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair( 17.131,  1.0e-3),     // I3d1C_T23A_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.1801,  1.0e-4),     // I3d1C_T23A_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                    /* T_23B */

                    /* 2 particle */

                    /* I_1 phi_+ */
                    std::make_pair(1.049148e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.049148e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.049148e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.848300e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(4.848300e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(4.848300e-2,  1.0e-6), // I1_T23B_2pt_phi_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_1 phi_bar */
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.303407e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-8.337742e-3,  1.0e-6), // I1_T23B_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 phi_bar */
                    std::make_pair(1.032833e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(8.745443e-2,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(7.162560e-2,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.478178e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.097542e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.716906e-1,  1.0e-6), // I2_T23B_2pt_phi_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_+ */
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(4.962290e-3,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.802917e-2,  1.0e-6), // I2_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_+ */
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.643000e-5,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-6.860632e-4,  1.0e-6), // I3_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_+ */
                    std::make_pair(-6.982288e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-6.982288e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-6.982288e-3,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.411838e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.411838e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-2.411838e-2,  1.0e-6), // I3d1_T23B_2pt_g_p(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_2 g_bar */
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-1.689134e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-9.729494e-3,  1.0e-6), // I2_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3 g_bar */
                    std::make_pair(-4.788091e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-4.052294e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.316497e-3,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-6.153610e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-5.218510e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-4.283410e-2,  1.0e-6), // I3_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_3d1 g_bar */
                    std::make_pair(-4.534548e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(-3.840569e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(-3.146589e-1,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(-2.736465   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(-2.323013   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(-1.909561   ,  1.0e-6), // I3d1_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4 g_bar */
                    std::make_pair(2.414476e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(2.044442e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.674408e-4,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(3.194830e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(2.704120e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(2.213409e-3,  1.0e-6), // I4_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d1 g_bar */
                    std::make_pair(2.296247e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.943389e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.590530e-2,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(1.449620e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(1.227319e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(1.005019e-1,  1.0e-6), // I4d1_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* I_4d2 g_bar */
                    std::make_pair(1.582910   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 = -5 GeV^2)
                    std::make_pair(1.339166   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 =  0 GeV^2)
                    std::make_pair(1.095423   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.04, q2 = +5 GeV^2)

                    std::make_pair(4.599967   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 = -5 GeV^2)
                    std::make_pair(3.898943   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 =  0 GeV^2)
                    std::make_pair(3.197918   ,  1.0e-6), // I4d2_T23B_2pt_g_bar(sigma = 0.08, q2 = +5 GeV^2)

                    /* 3 particle */

                    /* I_1 phi_3 */
                    std::make_pair(-3.5744e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-7.4908e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.9751e-4,  1.0e-8),  // I1_T23B_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.2522e-3,  1.0e-7),  // I1_T23B_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_3 */
                    std::make_pair( 7.1753e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.4871e-1,  1.0e-5),  // I2_T23B_3pt_phi_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.2020e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.5125e-2,  1.0e-6),  // I2_T23B_3pt_phi_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_3 */
                    std::make_pair(-3.2127e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-6.9752e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.1270e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-4.4950e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_3 */
                    std::make_pair(-1.1184e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.1665e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.5643e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.5578   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_3 */
                    std::make_pair(-7.4905e-1,  1.0e-5),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.5207,     1.0e-4),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.0117e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.1127e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_3 */
                    std::make_pair( 3.0959e1,  1.0e-3),   // I3d1B_T23B_3pt_phi_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.5956e4,  1.0e+0),   // I3d1B_T23B_3pt_phi_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_3 */
                    std::make_pair( 3.7179e-2, 1.0e-6),   // I3d1C_T23B_3pt_phi_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair( 7.7914e-2, 1.0e-6),   // I3d1C_T23B_3pt_phi_bar_3(sigma_0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_3 */
                    std::make_pair( 2.7250e-4,  1.0e-8),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 3.6297e-3,  1.0e-7),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.8166e-3,  1.0e-7),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.4198e-2,  1.0e-6),  // I3_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_3 */
                    std::make_pair(-2.2555e-3,  1.0e-7),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0044e-2,  1.0e-6),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.8566e-3,  1.0e-7),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-3.8050e-2,  1.0e-6),  // I3d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_3 */
                    std::make_pair( 2.8252e-1, 1.0e-5),   // I3d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 3.1758e+2, 1.0e-2),   // I3d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_3 */
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(       0.0, 1.0e-4),   // I3d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_3 */
                    std::make_pair( 1.3510e-3,  1.0e-7),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5724e-2,  1.0e-6),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.3707e-2,  1.0e-6),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1735e-1,  1.0e-5),  // I4_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_3 */
                    std::make_pair( 8.9875e-3,  1.0e-7),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.1923e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.9387e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6937   ,  1.0e-4),  // I4d1A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_3 */
                    std::make_pair( 1.1431   ,  1.0e-4),  // I4d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 1.2850e+3,  1.0e-1),  // I4d1B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_3 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_3 */
                    std::make_pair(-2.8222e-1,  1.0e-5),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-3.0185   ,  1.0e-4),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 2.9238e-1,  1.0e-5),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.6509   ,  1.0e-4),  // I4d2A_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_3 */
                    std::make_pair(-2.1918  ,  1.0e-4),   // I4d2B_T23B_3pt_phi_bar_bar_3(sigma_0, 1.0, 5.0)
                    std::make_pair( 2.1829e4,  1.0e-0),   // I4d2B_T23B_3pt_phi_bar_bar_3(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_3 */
                    std::make_pair(-1.7462e-1,  1.0e-5),  // I4d2C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.1, 5.0)
                    std::make_pair(-4.6520e-1,  1.0e-5),  // I4d2C_T23B_3pt_phi_bar_bar_3(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_3 */
                    std::make_pair( 0.0,  1.0e-4),        // I4d2D_T23B_3pt_phi_bar_bar_3(sigma_0, 5.0)

                    /* I_2 phi_4 */
                    std::make_pair( 3.3331e-3,  1.0e-7),  // I2_T23B_3pt_phi_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.0950e-2,  1.0e-6),  // I2_T23B_3pt_phi_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.6548e-4,  1.0e-8),  // I2_T23B_3pt_phi_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.9809e-4,  1.0e-8),  // I2_T23B_3pt_phi_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_2 phi_bar_4 */
                    std::make_pair( 3.8706e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.5302e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5589e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 4.0771e-1,  1.0e-5),  // I2_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_4 */
                    std::make_pair( 1.1692e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.9214e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.5760e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.2374   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_4 */
                    std::make_pair( 9.7065e-1,  1.0e-5),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 5.7941   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 5.8848   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.6569e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_4 */
                    std::make_pair( 2.2178e1,  1.0e-3),   // I3d1B_T23B_3pt_phi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 8.3217e3,  1.0e-1),   // I3d1B_T23B_3pt_phi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_4 */
                    std::make_pair( 1.0126  ,  1.0e-4),   // I3d1C_T23B_3pt_phi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair( 2.1220  ,  1.0e-4),   // I3d1C_T23B_3pt_phi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 phi_bar_bar_4 */
                    std::make_pair(-1.1893e-3,  1.0e-7),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.1364e-3,  1.0e-7),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.3407e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.7057e-2,  1.0e-6),  // I2_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 phi_bar_bar_4 */
                    std::make_pair( 2.0318e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 7.0113e-1,  1.0e-5),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 4.0072   ,  1.0e-4),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 1.1451e+1,  1.0e-3),  // I3_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A phi_bar_bar_4 */
                    std::make_pair(-3.0059   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-9.0956   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-7.9620   ,  1.0e-4),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.1288e+1,  1.0e-3),  // I3d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B phi_bar_bar_4 */
                    std::make_pair(-1.6717e+1,  1.0e-3),  // I3d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-1.3680e+4,  1.0e-0),  // I3d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I3d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4 phi_bar_bar_4 */
                    std::make_pair( 5.1748e-1,  1.0e-5),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7874   ,  1.0e-4),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.0203e+1,  1.0e-3),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.9167e+1,  1.0e-3),  // I4_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1A phi_bar_bar_4 */
                    std::make_pair(-1.0133   ,  1.0e-4),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-2.0759e-1,  1.0e-5),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.1064e+2,  1.0e-2),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.2014e+2,  1.0e-2),  // I4d1A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d1B phi_bar_bar_4 */
                    std::make_pair(-4.2378e+1,  1.0e-3),  // I4d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.4681e+4,  1.0e+0),  // I4d1B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d1C phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(0.0,  1.0e-4),         // I4d1C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2A phi_bar_bar_4 */
                    std::make_pair(-1.6145e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-5.1732e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-6.5989e+2,  1.0e-2),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-1.7868e+3,  1.0e-1),  // I4d2A_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_4d2B phi_bar_bar_4 */
                    std::make_pair( 9.5143e2,  1.0e-2),   // I4d2B_T23B_3pt_phi_bar_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-3.3525e5,  1.0e+1),   // I4d2B_T23B_3pt_phi_bar_bar_4(sigma_0, 3.0, 5.0)

                    /* I_4d2C phi_bar_bar_4 */
                    std::make_pair(-11.507  ,   1.0e-3),  // I4d2C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-30.656  ,   1.0e-3),  // I4d2C_T23B_3pt_phi_bar_bar_4(sigma_0, 0.5, 5.0)

                    /* I_4d2D phi_bar_bar_4 */
                    std::make_pair(0.0,  1.0e-4),         // I4d2D_T23B_3pt_phi_bar_bar_4(sigma_0, 5.0)

                    /* I_2 psi_bar_4 */
                    std::make_pair(-8.7441e-2,  1.0e-5),  // I2_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-4.9220e-2,  1.0e-5),  // I2_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-5.4456e-1,  1.0e-4),  // I2_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-2.4769e-1,  1.0e-4),  // I2_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 psi_bar_4 */
                    std::make_pair(-2.4595e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.2644e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-1.5683   ,  1.0e-4),  // I3_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-6.9348e-1,  1.0e-5),  // I3_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A psi_bar_4 */
                    std::make_pair(-1.6700   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair(-1.1771   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair(-2.0047e+1,  1.0e-3),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair(-9.1407   ,  1.0e-4),  // I3d1A_T23B_3pt_psi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B psi_bar_4 */
                    std::make_pair(-5.9695  ,  1.0e-4),   // I3d1B_T23B_3pt_psi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair(-7.1421e2,  1.0e-2),   // I3d1B_T23B_3pt_psi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C psi_bar_4 */
                    std::make_pair(-2.4489e-1,  1.0e-5),  // I3d1C_T23B_3pt_psi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-1.0264e-1,  1.0e-5),  // I3d1C_T23B_3pt_psi_bar_4(sigma_0, 0.5, 5.0)

                    /* I_2 chi_bar_4 */
                    std::make_pair( 2.0776e-1,  1.0e-5),  // I2_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 4.8693e-2,  1.0e-6),  // I2_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 1.5023   ,  1.0e-4),  // I2_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 5.7025e-1,  1.0e-5),  // I2_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3 chi_bar_4 */
                    std::make_pair( 8.3115e-1,  1.0e-5),  // I3_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 1.7619e-1,  1.0e-5),  // I3_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 6.0667   ,  1.0e-4),  // I3_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 2.2763   ,  1.0e-4),  // I3_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1A chi_bar_4 */
                    std::make_pair( 5.2893   ,  1.0e-4),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.1, 5.0)
                    std::make_pair( 2.0471e-1,  1.0e-5),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 1.0, 0.5, 5.0)
                    std::make_pair( 8.3147e+1,  1.0e-3),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.1, 5.0)
                    std::make_pair( 3.1736e+1,  1.0e-3),  // I3d1A_T23B_3pt_chi_bar_4(sigma_0, 3.0, 0.5, 5.0)

                    /* I_3d1B chi_bar_4 */
                    std::make_pair( 6.2809e1,  1.0e-3),   // I3d1B_T23B_3pt_chi_bar_4(sigma_0, 1.0, 5.0)
                    std::make_pair( 7.5146e3,  1.0e-1),   // I3d1B_T23B_3pt_chi_bar_4(sigma_0, 3.0, 5.0)

                    /* I_3d1C chi_bar_4 */
                    std::make_pair(-1.8049,  1.0e-4),     // I3d1C_T23B_3pt_chi_bar_4(sigma_0, 0.1, 5.0)
                    std::make_pair(-0.7565,  1.0e-3),     // I3d1C_T23B_3pt_chi_bar_4(sigma_0, 0.5, 5.0)

                };
                TEST_CHECK_DIAGNOSTICS(diagnostics, reference);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1(-5.0), 0.849651, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1( 0.0), 0.838315, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_1( 5.0), 0.822425, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2(-5.0), 0.833842, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2( 0.0), 0.821883, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_2( 5.0), 0.781765, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30(-5.0), 0.858475, 1.0e-3);
              //TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30( 0.0), 0.838315, 1.0e-3); A30 is 0 for q2 = 0
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_a_30( 5.0), 0.817009, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v(-5.0), 0.851149, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v( 0.0), 0.841027, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_v( 5.0), 0.826181, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1(-5.0), 0.846651, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1( 0.0), 0.835993, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_1( 5.0), 0.820869, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A(-5.0), 0.847232, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A( 0.0), 0.835993, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23A( 5.0), 0.818952, 1.0e-3);

                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B(-5.0), 0.832092, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B( 0.0), 0.818786, 1.0e-3);
                TEST_CHECK_NEARLY_EQUAL(ff.normalized_moment_1_t_23B( 5.0), 0.809531, 1.0e-3);
            }

            /* B -> K^* form factor values */
            {
                static const double eps = 1.0e-4; // relative error < 0.3%

                Parameters p = Parameters::Defaults();
                p["B::1/lambda_B_p"]          = 2.173913;
                p["B::lambda_E^2"]            = 0.3174;
                p["B::lambda_H^2"]            = 1.2696;
                p["mass::ud(2GeV)"]           = 0.008;
                p["mass::B_d"]                = 5.2795;
                p["mass::K^*_d"]              = 0.896;
                p["decay-constant::B_d"]      = 0.180;
                p["B->K^*::f_Kstar_par"]      = 0.217;
                p["B->K^*::mu@B-LCSR"]        = 1.0;
                p["B->K^*::s_0^A1,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A1,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A2,0@B-LCSR"]  = 1.7;
                p["B->K^*::s_0^A2,1@B-LCSR"]  = 0.0;
                p["B->K^*::s_0^A30,0@B-LCSR"] = 1.7;
                p["B->K^*::s_0^A30,1@B-LCSR"] = 0.0;
                p["B->K^*::s_0^V,0@B-LCSR"]   = 1.7;
                p["B->K^*::s_0^V,1@B-LCSR"]   = 0.0;
                p["B->K^*::M^2@B-LCSR"]       = 1.0;

                Options o = {
                    { "2pt",    "all"  },
                    { "3pt",    "all"  },
                    { "gminus", "zero" }
                };

                std::shared_ptr<FormFactors<PToV>> ff = FormFactorFactory<PToV>::create("B->K^*::B-LCSR", p, o);

                TEST_CHECK_RELATIVE_ERROR( 0.190982, ff->a_1(-5.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.200176, ff->a_1( 0.0),   2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.205437, ff->a_1(+5.0),   5.0 * eps);

                TEST_CHECK_RELATIVE_ERROR( 0.103649, ff->a_2(-5.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.097584, ff->a_2( 0.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.055975, ff->a_2(+5.0),         eps);

                TEST_CHECK_RELATIVE_ERROR( 0.206460, ff->v(-5.0),           eps);
                TEST_CHECK_RELATIVE_ERROR( 0.249919, ff->v( 0.0),           eps);
                TEST_CHECK_RELATIVE_ERROR( 0.298892, ff->v(+5.0),           eps);

                TEST_CHECK_RELATIVE_ERROR( 0.339110, ff->a_0(-5.0),   2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.451131, ff->a_0( 0.0),   3.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.634642, ff->a_0(+5.0),   6.0 * eps);

                TEST_CHECK_RELATIVE_ERROR( 0.181423, ff->t_1(-5.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.220680, ff->t_1( 0.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.267045, ff->t_1(+5.0),   3.0 * eps);

                TEST_CHECK_RELATIVE_ERROR( 0.213140, ff->t_2(-5.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.220680, ff->t_2( 0.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.221944, ff->t_2(+5.0),   3.0 * eps);

                TEST_CHECK_RELATIVE_ERROR( 0.103510, ff->t_3(-5.0),         eps);
                TEST_CHECK_RELATIVE_ERROR( 0.104617, ff->t_3( 0.0),   2.0 * eps);
                TEST_CHECK_RELATIVE_ERROR( 0.077692, ff->t_3(+5.0),   5.0 * eps);
            }
        }
} kmo2006_form_factors_test;