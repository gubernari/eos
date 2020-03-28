/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Nico Gubernari
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
#include <eos/form-factors/form-factors.hh>
#include <eos/form-factors/baryonic-hqet.hh>
#include <eos/form-factors/baryonic-impl.hh>

#include <vector>

using namespace test;
using namespace eos;

class LbToLcHQETFormFactorsTest :
    public TestCase
{
    public:
        LbToLcHQETFormFactorsTest() :
            TestCase("variables_test")
        {
        }

        virtual void run() const
        {
            static const double eps = 1.0e-4;


            {
                Parameters p = Parameters::Defaults();
                p["mass::B_d"]                              =  5.27942;
                p["mass::D_d"]                              =  1.86723;
                p["mass::Lambda_b"]                         =  5.6194;
                p["mass::Lambda_c"]                         =  2.2865;
                p["Lambda_b->Lambda_c::zeta'(1)@HQET"]      =  1.2;
                p["Lambda_b->Lambda_c::zeta''(1)@HQET"]     =  2.3;
                p["Lambda_b->Lambda_c::zeta'''(1)@HQET"]    =  3.4;
                p["Lambda_b->Lambda_c::zeta''''(1)@HQET"]   =  4.5;
                p["Lambda_b->Lambda_c::zeta'''''(1)@HQET"]  =  5.6;

                auto oo = Options{
                    { "z-order-lp",   "5" },
                    { "z-order-slp",  "2" },
                    { "z-order-sslp", "1" }
                };

                HQETVariables ff(p, oo);//what should I put here?
                Diagnostics diag = ff.diagnostics();

                static const std::vector<std::pair<double, double>> ref
                {
                    /* Inputs */
                    std::make_pair(+11.1082, eps), // q2(w = 1.0)
                    std::make_pair(+51.0746, eps), // sp
                    std::make_pair(+0.0, eps), // z(q2(w = 1))
                    std::make_pair(+1.68446, eps), // zeta_power_series(q2(w = 1.0))
                    std::make_pair(+1.59472, eps), // zeta_power_series(q2(w = 2.0))
                    std::make_pair(+1.22334, eps), // zeta_power_series(q2(w = 7.0))
                };

                TEST_CHECK_DIAGNOSTICS(diag, ref);
            }
        }
} variables_test;
