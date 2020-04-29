/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2015, 2016, 2017 Danny van Dyk
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
#include <eos/observable.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/wilson-polynomial.hh>

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

// Uncomment the following #define to generate new test data for the Bobeth compatibility tests
//#define EOS_GENERATE_TEST_DATA

using namespace test;
using namespace eos;

class BToKstarDileptonLowRecoilTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_test")
        {
        }

        virtual void run() const
        {
            /* Low Recoil (SM) */
            {
                Parameters p = Parameters::Defaults();
                p["life_time::B_d"] = 1.530e-12;
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = -0.331;
                p["b->s::c8"] = -0.18100000;
                p["b->smumu::Re{c9}"] = +4.27;
                p["b->smumu::Re{c10}"] = -4.173;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // muon mass near zero to avoid artificial divergence
                p["mass::mu"] = 1e-5;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("tag", "GP2004");
                oo.set("l", "mu");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton d(p, oo);

                /* q^2 = [14.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(14.00, 19.21), -0.4093, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(14.00, 19.21),  +0.3497, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(14.00, 19.21),     -0.4835, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(14.00, 19.21),     +1.6893, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(14.00, 19.21),     +0.5758, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(14.00, 19.21),      0.1244, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(14.00, 19.21),    -0.8391, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(14.00, 19.21),     0.0,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(14.00, 19.21),                        +0.9967, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(14.00, 19.21),                        -0.9727, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(14.00, 19.21),                        -0.9587, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(14.00, 19.21),                         0.0,    eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(14.00, 19.21),                         0.0,    eps);
                }

                /* q^2 = [16.00, 19.21] */
                {
                    const double eps = 1e-4;

                    TEST_CHECK_NEARLY_EQUAL(d.integrated_forward_backward_asymmetry(16.00, 19.21), -0.381708, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_longitudinal_polarisation(16.00, 19.21),  +0.337697, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_2(16.00, 19.21),     -0.599389, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_3(16.00, 19.21),     +1.99535,  eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_4(16.00, 19.21),     +0.486256, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_5(16.00, 19.21),      0.112158, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_re(16.00, 19.21),    -0.768382, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_transverse_asymmetry_im(16.00, 19.21),     0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_1(16.00, 19.21),                        +0.998622, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_2(16.00, 19.21),                        -0.970214, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_3(16.00, 19.21),                        -0.959887, eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_4(16.00, 19.21),                         0.0,      eps);
                    TEST_CHECK_NEARLY_EQUAL(d.integrated_h_5(16.00, 19.21),                         0.0,      eps);
                }

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -9.860564941316e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), +8.071641897174e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), -2.941484608501e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +7.179697602811e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), -5.877142772730e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), +2.141760651448e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -1.139839686524e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -3.400232049605e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), +9.330477335285e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), -3.400232049605e-12, eps);
                }
            }

            /* Low Recoil (Benchmark Point) */
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::Im{c7}"] = -0.3;
                p["b->s::c8"] = -0.181;
                p["b->smumu::Re{c9}"] = 0.0;
                p["b->smumu::Im{c9}"] = 4.2;
                p["b->smumu::Re{c10}"] = 0.0;
                p["b->smumu::Im{c10}"] = -4.2;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;
                // mu mass
                p["mass::mu"] = 1e-5;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("tag", "GP2004");
                oo.set("l", "mu");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -2.41522826885e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -9.96805582174e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), -2.41522826886e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), +7.68695280669e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +7.25796411402e-11, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), +1.75858165484e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), -5.59704205262e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -1.15226517859e-10, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), -2.79190193386e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), +8.88579298412e-12, eps);
                }
            }

            /* Low Recoil (Zero Point for C_7 = C_9 = C_10 = 0) */
            {
                Parameters p = Parameters::Defaults();
                p["b->s::c1"] = -0.32300000;
                p["b->s::c2"] = +1.00931000;
                p["b->s::c3"] = -0.00522869;
                p["b->s::c4"] = -0.08794730;
                p["b->s::c5"] = +0.00037476;
                p["b->s::c6"] = +0.00105859;
                p["b->s::Re{c7}"] = 0.0;
                p["b->s::c8"] = -0.181;
                p["b->smumu::Re{c9}"] = 0.0;
                p["b->smumu::Re{c10}"] = 0.0;
                // PDG 2008 CKM parameters
                p["CKM::A"] = 0.814;
                p["CKM::lambda"] = 0.2257;
                p["CKM::rhobar"] = 0.135;
                p["CKM::etabar"] = 0.349;
                // Kaon mass
                p["mass::K^*_d"] = 0.896;
                // B mass
                p["mass::B_d"] = 5.27953;
                // s quark mass
                p["mass::s(2GeV)"] = 0.0;
                // b quark mass
                p["mass::b(MSbar)"] = 4.2;

                Options oo;
                oo.set("model", "WilsonScan");
                oo.set("tag", "GP2004");
                oo.set("l", "mu");
                oo.set("form-factors", "BZ2004");

                BToKstarDilepton d(p, oo);

                /* transversity amplitudes at q^2 = 16.00 GeV^2 */
                {
                    static const double eps = 1e-19; // 1e-7 smaller than results

                    auto amps = d.amplitudes(16.00);

                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_left),  -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_left),  -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_long_right), -2.413541335202e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_long_right), -2.939430107299e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_left),  +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_left),  +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_perp_right), +1.757353360762e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_perp_right), +2.140264723229e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_left),  -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_left),  -3.397857132935e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(real(amps.a_para_right), -2.789951909754e-12, eps);
                    TEST_CHECK_NEARLY_EQUAL(imag(amps.a_para_right), -3.397857132935e-12, eps);
                }
            }
        }
} b_to_kstar_dilepton_low_recoil_test;

class BToKstarDileptonLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKstarDileptonLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_kstar_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            static const std::vector<std::string> variation_names
            {
                "b->s::Re{c7}",      "b->s::Im{c7}",      "b->s::Re{c7'}",      "b->s::Im{c7'}",
                "b->smumu::Re{c9}",  "b->smumu::Im{c9}",  "b->smumu::Re{c9'}",  "b->smumu::Im{c9'}",
                "b->smumu::Re{c10}", "b->smumu::Im{c10}", "b->smumu::Re{c10'}", "b->smumu::Im{c10'}",
            };

            Parameters p = Parameters::Defaults();
            // comparison done for zero lepton mass
            // but this leads to a NaN in the timelike transversity amplitude
            // so make the mass very small
            p["mass::mu"] = 1e-5;
            Options o;
            o.set("model", "WilsonScan");
            o.set("l", "mu");
            o.set("form-factors", "BZ2004");

            std::vector<Parameter> variations;

            for (auto n = variation_names.cbegin(), n_end = variation_names.cend() ; n != n_end ; ++n)
            {
                variations.push_back(p[*n]);
            }

            Kinematics k;
            k.declare("s_min"); k.set("s_min", 14.18);
            k.declare("s_max"); k.set("s_max", 19.21);

            std::vector<ObservablePtr> observables;
            observables.push_back(Observable::make("B->K^*ll::BR,tag=GP2004,q=d,l=mu",   p, k, o));
            observables.push_back(Observable::make("B->K^*ll::A_FB,tag=GP2004,q=d,l=mu", p, k, o));
            observables.push_back(Observable::make("B->K^*ll::F_L,tag=GP2004,q=d,l=mu",  p, k, o));

            std::string filename(EOS_SRCDIR "/eos/rare-b-decays/exclusive-b-to-s-dilepton-low-recoil_TEST-btokstarll.data");
#ifdef EOS_GENERATE_TEST_DATA
            {
                std::cout << "-- GENERATING test case data for B->K^*ll@LowRecoil --" << std::endl;
                RandomNumberGenerator rng;
                std::fstream file(filename.c_str(), std::fstream::out);
                file.precision(17);

                for (int i = 0 ; i < 1000 ; ++i)
                {
                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        *v = v->min() + (v->max() - v->min()) * rng();

                        file << *v << '\t';
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        file << (*o)->evaluate() << '\t';
                    }
                    file << std::endl;
                }
            }
#else
            // Verify the test case data
            {
                std::cout << "-- Verifying test case data for B->K^*ll@LowRecoil --" << std::endl;
                std::fstream file(filename.c_str(), std::fstream::in);

                std::string line;
                while (file)
                {
                    std::getline(file, line);
                    if (line.empty())
                        break;

                    std::stringstream ss(line);

                    for (auto v = variations.begin(), v_end = variations.end() ; v != v_end ; ++v)
                    {
                        double value;
                        ss >> value;
                        *v = value;
                    }

                    for (auto o = observables.cbegin(), o_end = observables.cend() ; o != o_end ; ++o)
                    {
                        double reference;
                        ss >> reference;

                        TEST_CHECK_RELATIVE_ERROR(reference, (*o)->evaluate(), 1e-3);
                    }
                }
            }
#endif
        }
} b_to_kstar_dilepton_low_recoil_bobeth_compatibility_test;

class BToKstarDileptonTensorLowRecoilBobethCompatibilityTest :
    public TestCase
{
    public:
        BToKstarDileptonTensorLowRecoilBobethCompatibilityTest() :
            TestCase("b_to_kstar_tensor_dilepton_low_recoil_bobeth_compatibility_test")
        {
        }

        virtual void run() const
        {
            // Christoph uses \Delta C instead of C for C9, C10
            // important to agree to alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["b->s::c1"] = -0.3231323312;
            p["b->s::c2"] = 1.009301831;
            p["b->s::c3"] = -0.005233499106;
            p["b->s::c4"] = -0.08829686414;
            p["b->s::c5"] = 0.0003601965805;
            p["b->s::c6"] = 0.001020749573;
            p["b->s::Re{c7}"] = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"] = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"] = -0.1827530948;
            p["b->smumu::Re{c9}"] = 4.294489364 + 1;
            p["b->smumu::Im{c9}"] = 0.5;
            p["b->smumu::Re{c9'}"] = 2;
            p["b->smumu::Im{c9'}"] = 1.5;
            p["b->smumu::Re{c10}"] = -4.196294696 + 3;
            p["b->smumu::Im{c10}"] = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["b->smumu::Re{cS}"] = 0.5;
            p["b->smumu::Im{cS}"] = 1;
            p["b->smumu::Re{cS'}"] = 0.6;
            p["b->smumu::Im{cS'}"] = 1.1;
            p["b->smumu::Re{cP}"] = 0.7;
            p["b->smumu::Im{cP}"] = 1.2;
            p["b->smumu::Re{cP'}"] = 0.8;
            p["b->smumu::Im{cP'}"] = 1.3;
            p["b->smumu::Re{cT}"] = 0.9;
            p["b->smumu::Im{cT}"] = 1.4;
            p["b->smumu::Re{cT5}"] = -1.0;
            p["b->smumu::Im{cT5}"] = -1.5;

            p["mass::s(2GeV)"] = 0.095;

            // increase sensitivity to m_l^2/q^2 terms
            p["mass::mu"] = 1.5;

            Options oo;
            oo.set("model", "WilsonScan");
            oo.set("scan-mode", "cartesian");
            oo.set("tag", "GP2004");
            oo.set("form-factors", "KMPW2010");
            oo.set("l", "mu");
            oo.set("q", "d");

            static const double q2 = 14.0;

            {
                double eps = 7e-3;

                BToKstarDilepton d(p, oo);
                auto amps = d.amplitudes(q2);

                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.021407965e-11, -1.564297789e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.645626526e-11,  8.331646455e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.455049449e-11,  4.565517978e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.565928422e-11,  2.577481906e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.541145186e-11,  4.618243535e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.607110071e-11,  2.607248335e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-4;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.677697256e-10, -3.507403558e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.767698228e-12,  2.767698228e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.38060e-11,      3.70316e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.64511e-11,      3.96767e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.46932e-11,      2.28561e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.63258e-11,      2.44887e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.12340e-11,      4.6851e-11 ),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.81106e-11,      4.37276e-11),     eps);
            }

            {
                oo.set("cp-conjugate", "true");
                BToKstarDilepton d(p, oo);

                auto amps = d.amplitudes(q2);

                double eps = 7e-3;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.021407965e-11,  1.843164004e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 8.645626526e-11, -8.05278024e-11 ), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-7.455049449e-11, -8.452349138e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 3.565928422e-11, -2.966165022e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-7.541145186e-11, -8.549962337e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 3.607110071e-11, -3.000420215e-11), eps);

                // nearly identically implemented, only difference from alpha_s
                eps = 1e-4;
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,       complex<double>(-1.677697256e-10,  3.507403558e-11), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,       complex<double>( 2.767698228e-12, -2.767698228e-12), eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp,  complex<double>( 2.3806e-11,      -3.70316e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long,  complex<double>( 2.64511e-11,     -3.96767e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp,  complex<double>( 1.46932e-11,     -2.28561e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp,  complex<double>( 1.63258e-11,     -2.44887e-11),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para,  complex<double>( 3.1234e-11,      -4.6851e-11 ),     eps);
                TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para,  complex<double>( 2.81106e-11,     -4.37276e-11),     eps);
            }
        }
} b_to_kstar_dilepton_tensor_low_recoil_bobeth_compatibility_test;
