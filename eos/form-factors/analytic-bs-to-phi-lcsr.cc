/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2018 Danny van Dyk
 * Copyright (c) 2019 Nico Gubernari
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

#include <eos/form-factors/analytic-b-to-v-lcsr-impl.hh>

namespace eos
{
    namespace lcsr
    {
        struct BsToPhi
        {
            constexpr static const char * B      = "B_s";
            constexpr static const char * m_B    = "mass::B_s";
            constexpr static const char * f_B    = "decay-constant::B_s";
            constexpr static const char * V      = "phi";
            constexpr static const char * m_V    = "mass::phi";
            constexpr static const char * f_V    = "decay-constant::phi";
            constexpr static const char   q_v    = 's';
            constexpr static const char   q_s    = 's';
            constexpr static const double chi2   = 1.0;
        };

        // B -> phi
        constexpr const char * BsToPhi::B;
        constexpr const char * BsToPhi::m_B;
        constexpr const char * BsToPhi::f_B;
        constexpr const char * BsToPhi::V;
        constexpr const char * BsToPhi::m_V;
        constexpr const char * BsToPhi::f_V;
        constexpr const char   BsToPhi::q_s;
    }

    template class AnalyticFormFactorBToVLCSR<lcsr::BsToPhi>;
}