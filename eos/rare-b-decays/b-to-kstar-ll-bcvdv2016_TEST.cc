/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017 Danny van Dyk
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

#include <array>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

using namespace test;
using namespace eos;

class BToKstarDileptonBCvDV2016MauriCompatibilityTest :
    public TestCase
{
    public:
    BToKstarDileptonBCvDV2016MauriCompatibilityTest() :
        TestCase("b_to_kstar_dilepton_BCvDV2016_mauri_compatibility_test")
    {
    }

    virtual void run() const
    {
        //todo we need tests for large-recoil-ff=ABBBSW2008
        {
            Parameters p = Parameters::Defaults();

            p["B->K^*::alpha^A0_0@BSZ2015" ] = 0.369196;
            p["B->K^*::alpha^A0_1@BSZ2015" ] =-1.365840;
            p["B->K^*::alpha^A0_2@BSZ2015" ] = 0.128191;
            p["B->K^*::alpha^A1_0@BSZ2015" ] = 0.297250;
            p["B->K^*::alpha^A1_1@BSZ2015" ] = 0.392378;
            p["B->K^*::alpha^A1_2@BSZ2015" ] = 1.189160;
            p["B->K^*::alpha^A12_1@BSZ2015"] = 0.533638;
            p["B->K^*::alpha^A12_2@BSZ2015"] = 0.483166;
            p["B->K^*::alpha^V_0@BSZ2015"  ] = 0.376313;
            p["B->K^*::alpha^V_1@BSZ2015"  ] =-1.165970;
            p["B->K^*::alpha^V_2@BSZ2015"  ] = 2.424430;
            p["B->K^*::alpha^T1_0@BSZ2015" ] = 0.312055;
            p["B->K^*::alpha^T1_1@BSZ2015" ] =-1.008930;
            p["B->K^*::alpha^T1_2@BSZ2015" ] = 1.527200;
            p["B->K^*::alpha^T2_1@BSZ2015" ] = 0.496846;
            p["B->K^*::alpha^T2_2@BSZ2015" ] = 1.614310;
            p["B->K^*::alpha^T23_0@BSZ2015"] = 0.667412;
            p["B->K^*::alpha^T23_1@BSZ2015"] = 1.318120;
            p["B->K^*::alpha^T23_2@BSZ2015"] = 3.823340;

            p["B->K^*ccbar::Re{alpha_0^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_0^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_0^long}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_0^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_0^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_0^long}@BCvDV2016"] = 0.0;

            p["B->K^*ccbar::Re{alpha_1^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_1^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_1^long}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_1^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_1^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_1^long}@BCvDV2016"] = 0.0;

            p["B->K^*ccbar::Re{alpha_2^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_2^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Re{alpha_2^long}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_2^perp}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_2^para}@BCvDV2016"] = 0.0;
            p["B->K^*ccbar::Im{alpha_2^long}@BCvDV2016"] = 0.0;

            p["b->s::Re{c7}"]       = -0.33726;
            p["b->s::Im{c7}"]       =  0.00000;
            p["b->s::Re{c7'}"]      =  0.00000;
            p["b->s::Im{c7'}"]      =  0.00000;
            p["b->smumu::Re{c9}"]   =  4.2734;
            p["b->smumu::Im{c9}"]   =  0.0000;
            p["b->smumu::Re{c9'}"]  =  0.0000;
            p["b->smumu::Im{c9'}"]  =  0.0000;
            p["b->smumu::Re{c10}"]  = -4.1661;
            p["b->smumu::Im{c10}"]  =  0.0000;
            p["b->smumu::Re{c10'}"] =  0.0000;
            p["b->smumu::Im{c10'}"] =  0.0000;

            Options oo;
            oo.set("model",          "WilsonScan");
            oo.set("tag",            "BCvDV2016");
            oo.set("correlator",     "BCvDV2016ModelA");
            oo.set("form-factors",   "BSZ2015");
            oo.set("l",              "mu");
            oo.set("q",              "d");

            static const double eps = 1e-3;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            std::cout << "amps.a_long_left  = " << amps.a_long_left  << std::endl;
            std::cout << "amps.a_long_right = " << amps.a_long_right  << std::endl;
            std::cout << "amps.a_para_left  = " << amps.a_para_left  << std::endl;
            std::cout << "amps.a_para_right = " << amps.a_para_right  << std::endl;
            std::cout << "amps.a_perp_left  = " << amps.a_perp_left  << std::endl;
            std::cout << "amps.a_perp_right = " << amps.a_perp_right  << std::endl;

//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.128830563e-10, +5.066150805e-12), eps);
//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>(+3.997103415e-11, +3.448212880e-11), eps);
//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-3.313040371e-11, +1.695431809e-11), eps);
//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>(+5.932084339e-11, +3.474608174e-11), eps);
//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>(+4.322089783e-11, +2.377670837e-12), eps);
//            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>(+9.294873950e-11, +1.087965051e-10), eps);
       }
    }
} b_to_kstar_dilepton_BCvDV2016_mauri_compatibility_test;
