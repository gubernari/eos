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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_CHARMONIUM_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_CHARMONIUM_HH 1

#include <eos/utils/complex.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/utils/private_implementation_pattern.hh>

namespace eos
{
    /**
     * Decay: B -> K^* psi
     * with psi a narrow charmonium: psi = { J/psi, psi(2S) }
     */
    class BToKstarCharmonium :
        public ParameterUser,
        public PrivateImplementationPattern<BToKstarCharmonium>
    {
        public:
            ///@name Basic operations.
            ///@{

            BToKstarCharmonium(const Parameters & parameters, const Options & options);
            ~BToKstarCharmonium();

            ///@}

            ///@name Observables
            ///@{

            /// Branching ratio
            double branching_ratio() const;

            /// Polarzation fractions and relative phases of the amplitudes
            double perp_polarization() const;
            double para_polarization() const;
            double long_polarization() const;

            double delta_perp_long() const;
            double delta_para_long() const;

            /// Angular observables as detected in the decay B -> K^* psi (-> l^+ l^-)
            double S_1s_LHCb() const;
            double S_1c_LHCb() const;
            double S_3_LHCb() const;
            double S_4_LHCb() const;
            double S_8_LHCb() const;
            double S_9_LHCb() const;

            /// Auxilliary observables used to reduce the number of fit parameters
            double re_ratio_perp() const;
            double re_ratio_para() const;
            double re_ratio_long() const;
            double im_ratio_perp() const;
            double im_ratio_para() const;
            double im_ratio_long() const;
            double abs_ratio_perp() const;
            double abs_ratio_para() const;
            double abs_ratio_long() const;

            ///@}
    };
}

#endif
