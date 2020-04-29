/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2016, 2017 Danny van Dyk
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

class BToKstarDileptonABBBSW2008BobethCompatibilityTest :
    public TestCase
{
    public:
    BToKstarDileptonABBBSW2008BobethCompatibilityTest() :
        TestCase("b_to_kstar_dilepton_ABBBSW2008_bobeth_compatibility_test")
    {
    }

    virtual void run() const
    {
        //todo we need tests for large-recoil-ff=ABBBSW2008
        {
            Parameters p = Parameters::Defaults();
            p["b->s::c1"]      = -0.3231323312;
            p["b->s::c2"]      = 1.009301831;
            p["b->s::c3"]      = -0.005233499106;
            p["b->s::c4"]      = -0.08829686414;
            p["b->s::c5"]      = 0.0003601965805;
            p["b->s::c6"]      = 0.001020749573;
            p["b->s::Re{c7}"]  = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"]  = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"]      = -0.1827530948;
            p["b->smumu::Re{c9}"]   = 4.294489364 + 1;
            p["b->smumu::Im{c9}"]   = 0.5;
            p["b->smumu::Re{c9'}"]  = 2;
            p["b->smumu::Im{c9'}"]  = 1.5;
            p["b->smumu::Re{c10}"]  = -4.196294696 + 3;
            p["b->smumu::Im{c10}"]  = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;

            Options oo;
            oo.set("model",          "WilsonScan");
            oo.set("tag",            "ABBBSW2008");
            oo.set("qcdf-integrals", "mixed");
            oo.set("form-factors",   "KMPW2010");
            oo.set("l",              "mu");
            oo.set("q",              "d");

            static const double eps = 1e-3;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.128830563e-10, +5.066150805e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>(+3.997103415e-11, +3.448212880e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-3.313040371e-11, +1.695431809e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>(+5.932084339e-11, +3.474608174e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>(+4.322089783e-11, +2.377670837e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>(+9.294873950e-11, +1.087965051e-10), eps);
       }

       // scalar and tensor
       {
            // important to agree on alpha_s, can change values by 1%
            Parameters p = Parameters::Defaults();
            p["b->s::c1"]      = -0.3231323312;
            p["b->s::c2"]      = 1.009301831;
            p["b->s::c3"]      = -0.005233499106;
            p["b->s::c4"]      = -0.08829686414;
            p["b->s::c5"]      = 0.0003601965805;
            p["b->s::c6"]      = 0.001020749573;
            p["b->s::Re{c7}"]  = -0.3370422989 + 0.1;
            p["b->s::Im{c7}"]  = 0.2;
            p["b->s::Re{c7'}"] = 0.3;
            p["b->s::Im{c7'}"] = 0.4;
            p["b->s::c8"]      = -0.1827530948;
            p["b->smumu::Re{c9}"]   = 4.294489364 + 1;
            p["b->smumu::Im{c9}"]   = 0.5;
            p["b->smumu::Re{c9'}"]  = 2;
            p["b->smumu::Im{c9'}"]  = 1.5;
            p["b->smumu::Re{c10}"]  = -4.196294696 + 3;
            p["b->smumu::Im{c10}"]  = 2.5;
            p["b->smumu::Re{c10'}"] = 4;
            p["b->smumu::Im{c10'}"] = 3.5;
            p["b->smumu::Re{cS}"]   = 0.5;
            p["b->smumu::Im{cS}"]   = 1;
            p["b->smumu::Re{cS'}"]  = 0.6;
            p["b->smumu::Im{cS'}"]  = 1.1;
            p["b->smumu::Re{cP}"]   = 0.7;
            p["b->smumu::Im{cP}"]   = 1.2;
            p["b->smumu::Re{cP'}"]  = 0.8;
            p["b->smumu::Im{cP'}"]  = 1.3;
            p["b->smumu::Re{cT}"]   = 0.9;
            p["b->smumu::Im{cT}"]   = 1.4;
            p["b->smumu::Re{cT5}"]  = 1.0;
            p["b->smumu::Im{cT5}"]  = 1.5;

            p["mass::s(2GeV)"] = 0.12;

            Options oo;
            oo.set("model",          "WilsonScan");
            oo.set("tag",            "ABBBSW2008");
            oo.set("qcdf-integrals", "mixed");
            oo.set("form-factors",   "KMPW2010");
            oo.set("l",              "mu");
            oo.set("q",              "d");

            static const double eps = 2e-4;
            static const double q2 = 6.0;

            BToKstarDilepton d(p, oo);
            auto amps = d.amplitudes(q2);

            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_left,  complex<double>(-1.121022032e-10, 5.991646324e-12), eps);
            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_right, complex<double>( 4.333216107e-11, 3.590418466e-11), eps);
            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_left,  complex<double>( 4.353425835e-11, 5.287884397e-12), eps);
            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_perp_right, complex<double>( 9.326263994e-11, 1.117078741e-10), eps);
            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_left,  complex<double>(-4.176304959e-11, 1.651237347e-11), eps);
            //TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_right, complex<double>( 5.964843896e-11, 3.602848917e-11), eps);

            // nearly identically implemented, only difference from alpha_s
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time,      complex<double>(-2.247078271e-10, -6.370327589e-11), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_scal,      complex<double>( 2.185643583e-12,  2.185643583e-12), eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_para_perp, complex<double>( 1.47786e-11,      2.2989e-11),      eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_long, complex<double>(-1.64207e-11,     -2.46311e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_perp, complex<double>( 2.4322e-11,       3.78342e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_perp, complex<double>(-2.70244e-11,     -4.05366e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_time_para, complex<double>(-3.24769e-11,     -4.87154e-11),     eps);
            TEST_CHECK_RELATIVE_ERROR_C(amps.a_long_para, complex<double>( 2.92292e-11,      4.54677e-11),     eps);
       }
    }
} b_to_kstar_dilepton_ABBBSW2008_bobeth_compatibility_test;
