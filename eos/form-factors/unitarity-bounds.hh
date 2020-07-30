/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    class BGLCoefficients :
        public virtual ParameterUser,
        public PrivateImplementationPattern<BGLCoefficients>
    {
        public:
            BGLCoefficients(const Parameters &, const Options &);
            ~BGLCoefficients();

            // B -> D form factors (q=u,d)
            // {{{
            double V1_a0() const;
            double V1_a1() const;
            double V1_a2() const;

            double S1_a0() const;
            double S1_a1() const;
            double S1_a2() const;
            // }}}

            // B -> D^* form factor (q=u,d)
            // {{{
            double A1_a0() const;
            double A1_a1() const;
            double A1_a2() const;

            double A5_a0() const;
            double A5_a1() const;
            double A5_a2() const;

            double V4_a0() const;
            double V4_a1() const;
            double V4_a2() const;

            double P1_a0() const;
            double P1_a1() const;
            double P1_a2() const;
            // }}}

            // B^* -> D form factor (q=u,d)
            // {{{
            double P2_a0() const;
            double P2_a1() const;
            double P2_a2() const;

            double V5_a0() const;
            double V5_a1() const;
            double V5_a2() const;

            double A2_a0() const;
            double A2_a1() const;
            double A2_a2() const;

            double A6_a0() const;
            double A6_a1() const;
            double A6_a2() const;
            // }}}

            // B^* -> D^* form factor (q=u,d)
            // {{{
            double S2_a0() const;
            double S2_a1() const;
            double S2_a2() const;

            double S3_a0() const;
            double S3_a1() const;
            double S3_a2() const;

            double P3_a0() const;
            double P3_a1() const;
            double P3_a2() const;

            double V2_a0() const;
            double V2_a1() const;
            double V2_a2() const;

            double V3_a0() const;
            double V3_a1() const;
            double V3_a2() const;

            double V6_a0() const;
            double V6_a1() const;
            double V6_a2() const;

            double V7_a0() const;
            double V7_a1() const;
            double V7_a2() const;

            double A3_a0() const;
            double A3_a1() const;
            double A3_a2() const;

            double A4_a0() const;
            double A4_a1() const;
            double A4_a2() const;

            double A7_a0() const;
            double A7_a1() const;
            double A7_a2() const;
            // }}}

            // B_s -> D_s form factors
            // {{{
            double V1s_a0() const;
            double V1s_a1() const;
            double V1s_a2() const;

            double S1s_a0() const;
            double S1s_a1() const;
            double S1s_a2() const;
            // }}}

            // B_s -> D_s^* form factor
            // {{{
            double A1s_a0() const;
            double A1s_a1() const;
            double A1s_a2() const;

            double A5s_a0() const;
            double A5s_a1() const;
            double A5s_a2() const;

            double V4s_a0() const;
            double V4s_a1() const;
            double V4s_a2() const;

            double P1s_a0() const;
            double P1s_a1() const;
            double P1s_a2() const;
            // }}}

            // B_s^* -> D_s form factor
            // {{{
            double P2s_a0() const;
            double P2s_a1() const;
            double P2s_a2() const;

            double V5s_a0() const;
            double V5s_a1() const;
            double V5s_a2() const;

            double A2s_a0() const;
            double A2s_a1() const;
            double A2s_a2() const;

            double A6s_a0() const;
            double A6s_a1() const;
            double A6s_a2() const;
            // }}}

            // B_s^* -> D_s^* form factor
            // {{{
            double S2s_a0() const;
            double S2s_a1() const;
            double S2s_a2() const;

            double S3s_a0() const;
            double S3s_a1() const;
            double S3s_a2() const;

            double P3s_a0() const;
            double P3s_a1() const;
            double P3s_a2() const;

            double V2s_a0() const;
            double V2s_a1() const;
            double V2s_a2() const;

            double V3s_a0() const;
            double V3s_a1() const;
            double V3s_a2() const;

            double V6s_a0() const;
            double V6s_a1() const;
            double V6s_a2() const;

            double V7s_a0() const;
            double V7s_a1() const;
            double V7s_a2() const;

            double A3s_a0() const;
            double A3s_a1() const;
            double A3s_a2() const;

            double A4s_a0() const;
            double A4s_a1() const;
            double A4s_a2() const;

            double A7s_a0() const;
            double A7s_a1() const;
            double A7s_a2() const;
            // }}}

            // CLN parameter
            double rho2s() const;
    };

    /* Unitarity bound implemented as discussed in [BJvD2019] */
    class HQETUnitarityBounds :
        public virtual ParameterUser,
        public PrivateImplementationPattern<HQETUnitarityBounds>
    {
        public:
            HQETUnitarityBounds(const Parameters &, const Options &);
            ~HQETUnitarityBounds();

            // unitarity bounds as priors
            double bound_0p_prior() const;

            double bound_0m_prior() const;

            double bound_1p_prior() const;

            double bound_1m_prior() const;

            // unitarity bounds as pseudo observables
            double bound_0p() const;

            double bound_0m() const;

            double bound_1p() const;

            double bound_1m() const;
    };
}
