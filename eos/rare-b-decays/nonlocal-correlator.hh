/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
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

#ifndef EOS_GUARD_EOS_RARE_B_DECAYS_NONLOCAL_CORRELATOR_HH
#define EOS_GUARD_EOS_RARE_B_DECAYS_NONLOCAL_CORRELATOR_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/utils/complex.hh>
#include <eos/utils/diagnostics.hh>
#include <eos/utils/options.hh>
#include <eos/utils/parameters.hh>
#include <eos/rare-b-decays/nonlocal-correlator-fwd.hh>

#include <memory>
#include <string>

namespace eos
{
    /**
     * Provides the hadronic matrix element of the non-local operator T{ cbar gamma^mu c(x), C_1 O_1 + C_2 O_2 }.
     * We decompose this matrix element as in [BCvDV:2017A], eq. (4).
     */
    template <typename Transition_>
    class NonlocalCorrelator;

    template <typename Transition_>
    using NonlocalCorrelatorPtr = std::shared_ptr<NonlocalCorrelator<Transition_>>;

    /**
     * Pseudoobservable in order to expose the nonlocal correlator, see NonlocalCorrelator.
     */
    template <typename Process_, typename Transition_>
    class NonlocalCorrelatorObservable;

    // P -> P

    template <>
    class NonlocalCorrelator<nc::PToP> :
        public ParameterUser
    {
        protected:
            ///@name Stubs to throw InternalError whenever an implementation without residues is called.
            ///@{

            complex<double> jpsi_residues_not_implemented() const;
            complex<double> psi2s_residues_not_implemented() const;
            complex<double> moments_not_implemented() const;

            ///@}

        public:
            ///@name Basic operations
            ///@{

            virtual ~NonlocalCorrelator();

            ///@}

            ///@name Evaluate the correlator at arbitrary q2 values.
            ///@{

            virtual complex<double> H_plus(const double & q2) const = 0;

            ///@}

            ///@name Evaluate the first normalized moment of the correlator.
            ///@{

            virtual complex<double> normalized_moment_A(const double & /*q2*/) const { return moments_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the correlator on the J/psi pole.
            ///@{

            virtual complex<double> H_plus_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the correlator on the psi(2S) pole.
            ///@{

            virtual complex<double> H_plus_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            ///@}

            /// Factory method.
            static NonlocalCorrelatorPtr<nc::PToP> make(const QualifiedName & name, const Parameters & p, const Options & o);

            ///@name Internal diagnostics for unit tests
            ///@{
            virtual Diagnostics diagnostics() const = 0;
            ///@}
    };

    template <typename Process_>
    class NonlocalCorrelatorObservable<Process_, nc::PToP> :
        public ParameterUser,
        public PrivateImplementationPattern<NonlocalCorrelatorObservable<Process_, nc::PToP>>
    {
        public:
            ///@name Basic operations
            ///@{

            NonlocalCorrelatorObservable(const Parameters &, const Options &);
            ~NonlocalCorrelatorObservable();

            ///@}

            ///@name Correlator normalized to form factor as observable
            ///@{

            double re_ratio_plus(const double & q2) const;

            double im_ratio_plus(const double & q2) const;

            ///@}

            ///@name Correlator as observable
            ///@{

            double re_H_plus(const double & q2) const;

            double im_H_plus(const double & q2) const;

            double abs_H_plus(const double & q2) const;

            double arg_H_plus(const double & q2) const;

            ///@}

            ///@name First moment of the correlator as observable
            ///@{

            double re_normalized_moment_A(const double & q2) const;

            double im_normalized_moment_A(const double & q2) const;

            ///@}
    };
    extern template class NonlocalCorrelatorObservable<nc::BToK, nc::PToP>;

    // P -> V

    template <>
    class NonlocalCorrelator<nc::PToV> :
        public ParameterUser
    {
        protected:
            ///@name Stubs to throw InternalError whenever an implementation without residues is called.
            ///@{

            complex<double> jpsi_residues_not_implemented() const;
            complex<double> psi2s_residues_not_implemented() const;
            complex<double> moments_not_implemented() const;

            ///@}

        public:
            ///@name Basic operations
            ///@{

            virtual ~NonlocalCorrelator();

            ///@}

            ///@name Evaluate the correlator at arbitrary q2 values.
            ///@{

            virtual complex<double> H_perp(const double & q2) const = 0;

            virtual complex<double> H_para(const double & q2) const = 0;

            virtual complex<double> H_long(const double & q2) const = 0;

            ///@}

            ///@name Evaluate the first normalized moment of the correlator.
            ///@{

            virtual complex<double> normalized_moment_V1(const double & /*q2*/) const { return moments_not_implemented(); };

            virtual complex<double> normalized_moment_V2(const double & /*q2*/) const { return moments_not_implemented(); };

            virtual complex<double> normalized_moment_V23(const double & /*q2*/) const { return moments_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the correlator on the J/psi pole.
            ///@{

            virtual complex<double> H_perp_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            virtual complex<double> H_para_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            virtual complex<double> H_long_residue_jpsi() const { return jpsi_residues_not_implemented(); };

            ///@}

            ///@name Evaluate the residue of the correlator on the psi(2S) pole.
            ///@{

            virtual complex<double> H_perp_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            virtual complex<double> H_para_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            virtual complex<double> H_long_residue_psi2s() const { return psi2s_residues_not_implemented(); };

            ///@}

            /// Factory method.
            static NonlocalCorrelatorPtr<nc::PToV> make(const QualifiedName & name, const Parameters & p, const Options & o);

            ///@name Internal diagnostics for unit tests
            ///@{
            virtual Diagnostics diagnostics() const = 0;
            ///@}
    };

    template <typename Process_>
    class NonlocalCorrelatorObservable<Process_, nc::PToV> :
        public ParameterUser,
        public PrivateImplementationPattern<NonlocalCorrelatorObservable<Process_, nc::PToV>>
    {
        public:
            ///@name Basic operations
            ///@{

            NonlocalCorrelatorObservable(const Parameters &, const Options &);
            ~NonlocalCorrelatorObservable();

            ///@}

            ///@name Correlator normalized to form factor as observable
            ///@{

            double re_ratio_perp(const double & q2) const;

            double im_ratio_perp(const double & q2) const;

            double re_ratio_para(const double & q2) const;

            double im_ratio_para(const double & q2) const;

            double re_ratio_zero(const double & q2) const;

            double im_ratio_zero(const double & q2) const;

            ///@}

            ///@name Correlator as observable
            ///@{

            double re_H_perp(const double & q2) const;

            double im_H_perp(const double & q2) const;

            double abs_H_perp(const double & q2) const;

            double arg_H_perp(const double & q2) const;

            double re_H_para(const double & q2) const;

            double im_H_para(const double & q2) const;

            double abs_H_para(const double & q2) const;

            double arg_H_para(const double & q2) const;

            double re_H_long(const double & q2) const;

            double im_H_long(const double & q2) const;

            double abs_H_long(const double & q2) const;

            double arg_H_long(const double & q2) const;

            ///@}

            ///@name First moment of the correlator as observable
            ///@{

            double re_normalized_moment_V1(const double & q2) const;

            double im_normalized_moment_V1(const double & q2) const;

            double re_normalized_moment_V2(const double & q2) const;

            double im_normalized_moment_V2(const double & q2) const;

            double re_normalized_moment_V23(const double & q2) const;

            double im_normalized_moment_V23(const double & q2) const;

            ///@}
    };
    extern template class NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>;
}

#endif
