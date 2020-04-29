/* vim: set sw=4 sts=4 et tw=150 foldmethod=marker : */

/*
 * Copyright (c) 2019 Danny van Dyk
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

#include <eos/observable-impl.hh>
#include <eos/rare-b-decays/decays.hh>
#include <eos/rare-b-decays/b-to-kstar-charmonium.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/rare-b-decays/b-to-kstar-ll.hh>
#include <eos/rare-b-decays/exclusive-b-to-dilepton.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-dilepton.hh>
#include <eos/rare-b-decays/inclusive-b-to-s-gamma.hh>
#include <eos/rare-b-decays/lambda-b-to-lambda-dilepton.hh>
#include <eos/rare-b-decays/nonlocal-correlator.hh>
#include <eos/utils/concrete_observable.hh>

namespace eos
{
    // B_q -> l^+ l^-
    // {{{
    ObservableGroup
    make_b_to_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour.)",
            {
                make_observable("B_q->ll::BR", R"(\mathcal{B}(B_q \to \ell^+\ell^-))",
                        &BToDilepton::branching_ratio_time_zero),

                make_observable("B_q->ll::BR@Untagged", R"(\left\langle\mathcal{B}(B_q \to \ell^+\ell^-)\right\rangle)",
                        &BToDilepton::branching_ratio_untagged_integrated),

                make_observable("B_q->ll::A_DeltaGamma", R"(\mathcal{A}_{\Delta\Gamma}(B_q \to \ell^+\ell^-))",
                        &BToDilepton::cp_asymmetry_del_gamma),

                make_observable("B_q->ll::S", R"(\mathcal{S}(B_q \to \ell^+\ell^-))",
                        &BToDilepton::cp_asymmetry_mixing_S),

                make_observable("B_q->ll::eff_lifetime", R"(\langle\tau\rangle(B_q \to \ell^+\ell^-))",
                        &BToDilepton::effective_lifetime),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

#if 0
    // B_q -> P l^+l^-
    // {{{
    ObservableGroup
    make_b_to_p_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to P \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour.)",
            {
                // B -> K ll, Large Recoil
                make_observable("B->Kll::d^2Gamma@LargeRecoil",
                        R"(d^2\mathcal{\Gamma(\bar{B}\to \bar{K}\ell^+\ell^-)}/(dq^2\, d\cos\theta_\ell))",
                        &BToKDilepton<LargeRecoil>::two_differential_decay_width,
                        std::make_tuple("s", "cos(theta_l)")),

                make_observable("B->Kll::dBR/ds@LargeRecoil",
                        R"(d\mathcal{B}(\bar{B}\to \bar{K}\ell^+\ell^-)/dq^2)",
                        &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->Kll::F_H(q2)@LargeRecoil",
                        R"(F_\text{H}(\bar{B}\to \bar{K}\ell^+\ell^-)(q^2))",
                        &BToKDilepton<LargeRecoil>::differential_flat_term,
                        std::make_tuple("q2")),

                make_observable("B->Kll::A_FB(q2)@LargeRecoil",
                        R"(A_\text{FB}(\bar{B}\to \bar{K}\ell^+\ell^-)(q^2))",
                        &BToKDilepton<LargeRecoil>::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable_ratio("B->Kll::R_K(q2)@LargeRecoil",
                        R"(R_K(q^2))",
                        &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "mu" } },
                        &BToKDilepton<LargeRecoil>::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "e" } }
                        ),

                make_observable("B->Kll::BR@LargeRecoil",
                        R"(\mathcal{B}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::BRavg@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::integrated_branching_ratio_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_CP@LargeRecoil",
                        R"(A_\text{CP}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        &BToKDilepton<LargeRecoil>::integrated_cp_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::Gamma@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::F_H@LargeRecoil",
                        R"(F_\text{H}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        &BToKDilepton<LargeRecoil>::integrated_flat_term,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::F_Havg@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::integrated_flat_term_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_FB@LargeRecoil",
                        R"(A_\text{FB}(\bar{B}\to \bar{K}\ell^+\ell^-))",
                        &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_FBavg@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable_ratio("B->Kll::R_K@LargeRecoil",
                        R"(R_K)",
                        &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "l", "mu" } },
                        &BToKDilepton<LargeRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max"),
                        Options{ { "l", "e" } }
                        ),

                make_observable("B->Kll::R_K@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::integrated_ratio_muons_electrons,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::a_l@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::a_l,
                        std::make_tuple("q2")),

                make_observable("B->Kll::b_l@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::b_l,
                        std::make_tuple("q2")),

                make_observable("B->Kll::c_l@LargeRecoil",
                        &BToKDilepton<LargeRecoil>::c_l,
                        std::make_tuple("q2")),

                // B -> K ll, Low Recoil
                make_observable("B->Kll::d^2Gamma@LowRecoil",
                        &BToKDilepton<LowRecoil>::two_differential_decay_width,
                        std::make_tuple("s", "cos(theta_l)")),

                make_observable("B->Kll::dBR/ds@LowRecoil",
                        &BToKDilepton<LowRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->Kll::F_H(q2)@LowRecoil",
                        &BToKDilepton<LowRecoil>::differential_flat_term,
                        std::make_tuple("q2")),

                make_observable("B->Kll::A_FB(q2)@LowRecoil",
                        &BToKDilepton<LowRecoil>::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B->Kll::R_K(q2)@LowRecoil",
                        &BToKDilepton<LowRecoil>::differential_ratio_muons_electrons,
                        std::make_tuple("q2")),

                make_observable("B->Kll::BR@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::BRavg@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_branching_ratio_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_CP@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_cp_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::Gamma@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::F_H@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_flat_term,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::F_Havg@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_flat_term_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_FB@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::A_FBavg@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_forward_backward_asymmetry_cp_averaged,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::R_K@LowRecoil",
                        &BToKDilepton<LowRecoil>::integrated_ratio_muons_electrons,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->Kll::a_l@LowRecoil",
                        &BToKDilepton<LowRecoil>::a_l,
                        std::make_tuple("q2")),

                make_observable("B->Kll::b_l@LowRecoil",
                        &BToKDilepton<LowRecoil>::b_l,
                        std::make_tuple("q2")),

                make_observable("B->Kll::c_l@LowRecoil",
                        &BToKDilepton<LowRecoil>::c_l,
                        std::make_tuple("q2")),

                make_observable("B->Kll::Re{c9eff}@LowRecoil",
                        &BToKDilepton<LowRecoil>::real_c9eff,
                        std::make_tuple("q2")),

                make_observable("B->Kll::Im{c9eff}@LowRecoil",
                        &BToKDilepton<LowRecoil>::imag_c9eff,
                        std::make_tuple("q2")),

                make_observable("B->Kll::Re{c7eff}@LowRecoil",
                        &BToKDilepton<LowRecoil>::real_c7eff,
                        std::make_tuple("q2")),

                make_observable("B->Kll::Im{c7eff}@LowRecoil",
                        &BToKDilepton<LowRecoil>::imag_c7eff,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}
#endif

    // B_q -> V gamma
    // {{{
    ObservableGroup
    make_b_to_v_gamma_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to V gamma$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour.)",
            {
                // B -> K^* gamma
                make_observable("B->K^*gamma::BR",
                        &BToKstarGamma::branching_ratio),

                make_observable("B->K^*gamma::BRavg",
                        &BToKstarGamma::branching_ratio_cp_averaged),

                make_observable("B->K^*gamma::A_CP",
                        &BToKstarGamma::cp_asymmetry),

                make_observable("B->K^*gamma::S_K^*gamma",
                        &BToKstarGamma::s_kstar_gamma),

                make_observable("B->K^*gamma::C_K^*gamma",
                        &BToKstarGamma::c_kstar_gamma),

                make_observable("B->K^*gamma::A_I",
                        &BToKstarGamma::isospin_asymmetry),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B_q -> V l^+l^-
    // {{{
    ObservableGroup
    make_b_to_v_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B_q \to V \ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour.)",
            {
                // B -> K^* ll, Large Recoil
                make_observable("B->K^*ll::d^4Gamma",
                        &BToKstarDilepton::decay_width,
                        std::make_tuple("s", "cos(theta_l)", "cos(theta_k)", "phi")),

                make_observable("B->K^*ll::dBR/ds",
                        &BToKstarDilepton::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_FB(q2)",
                        &BToKstarDilepton::differential_forward_backward_asymmetry,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^2(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^3(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^4(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^5(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^re(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_re,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::A_T^im(q2)",
                        &BToKstarDilepton::differential_transverse_asymmetry_im,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::P'_4(q2)",
                        &BToKstarDilepton::differential_p_prime_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::P'_5(q2)",
                        &BToKstarDilepton::differential_p_prime_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::P'_6(q2)",
                        &BToKstarDilepton::differential_p_prime_6,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::F_L(q2)",
                        &BToKstarDilepton::differential_longitudinal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::F_T(q2)",
                        &BToKstarDilepton::differential_transversal_polarisation,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_1s(q2)",
                        &BToKstarDilepton::differential_j_1s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_1c(q2)",
                        &BToKstarDilepton::differential_j_1c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_2s(q2)",
                        &BToKstarDilepton::differential_j_2s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_2c(q2)",
                        &BToKstarDilepton::differential_j_2c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_3(q2)",
                        &BToKstarDilepton::differential_j_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_4(q2)",
                        &BToKstarDilepton::differential_j_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_5(q2)",
                        &BToKstarDilepton::differential_j_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_6s(q2)",
                        &BToKstarDilepton::differential_j_6s,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_6c(q2)",
                        &BToKstarDilepton::differential_j_6c,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_7(q2)",
                        &BToKstarDilepton::differential_j_7,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_8(q2)",
                        &BToKstarDilepton::differential_j_8,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::J_9(q2)",
                        &BToKstarDilepton::differential_j_9,
                        std::make_tuple("q2")),

#if 0
                make_observable("B->K^*ll::D_4(q2)",
                        &BToKstarDilepton::differential_d_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::D_5(q2)",
                        &BToKstarDilepton::differential_d_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::D_6s(q2)",
                        &BToKstarDilepton::differential_d_6s,
                        std::make_tuple("q2")),
#endif

                make_observable_ratio("B->K^*ll::R_K^*(q2)",
                        R"(R_{K^*}(q^2))",
                        &BToKstarDilepton::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "mu" } },
                        &BToKstarDilepton::differential_branching_ratio,
                        std::make_tuple("q2"),
                        Options{ { "l", "e" } }
                        ),

                make_observable("B->K^*ll::A_FB",
                        &BToKstarDilepton::integrated_forward_backward_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::BR",
                        &BToKstarDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_CP",
                        &BToKstarDilepton::integrated_cp_asymmetry,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::F_L",
                        &BToKstarDilepton::integrated_longitudinal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::F_T",
                        &BToKstarDilepton::integrated_transversal_polarisation,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^2",
                        &BToKstarDilepton::integrated_transverse_asymmetry_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^3",
                        &BToKstarDilepton::integrated_transverse_asymmetry_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^4",
                        &BToKstarDilepton::integrated_transverse_asymmetry_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^5",
                        &BToKstarDilepton::integrated_transverse_asymmetry_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^re",
                        &BToKstarDilepton::integrated_transverse_asymmetry_re,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::A_T^im",
                        &BToKstarDilepton::integrated_transverse_asymmetry_im,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::P'_4",
                        &BToKstarDilepton::integrated_p_prime_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::P'_5",
                        &BToKstarDilepton::integrated_p_prime_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::P'_6",
                        &BToKstarDilepton::integrated_p_prime_6,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^1(q2)",
                        &BToKstarDilepton::differential_h_1,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^2(q2)",
                        &BToKstarDilepton::differential_h_2,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^3(q2)",
                        &BToKstarDilepton::differential_h_3,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^4(q2)",
                        &BToKstarDilepton::differential_h_4,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^5(q2)",
                        &BToKstarDilepton::differential_h_5,
                        std::make_tuple("q2")),

                make_observable("B->K^*ll::H_T^1",
                        &BToKstarDilepton::integrated_h_1,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^2",
                        &BToKstarDilepton::integrated_h_2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^3",
                        &BToKstarDilepton::integrated_h_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^4",
                        &BToKstarDilepton::integrated_h_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::H_T^5",
                        &BToKstarDilepton::integrated_h_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::s_0^A_FB",
                        &BToKstarDilepton::a_fb_zero_crossing),

                make_observable("B->K^*ll::Gamma",
                        &BToKstarDilepton::integrated_decay_width,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_1s",
                        &BToKstarDilepton::integrated_j_1s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_1c",
                        &BToKstarDilepton::integrated_j_1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_2s",
                        &BToKstarDilepton::integrated_j_2s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_2c",
                        &BToKstarDilepton::integrated_j_2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_3",
                        &BToKstarDilepton::integrated_j_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_4",
                        &BToKstarDilepton::integrated_j_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_5",
                        &BToKstarDilepton::integrated_j_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_6s",
                        &BToKstarDilepton::integrated_j_6s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_6c",
                        &BToKstarDilepton::integrated_j_6c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_7",
                        &BToKstarDilepton::integrated_j_7,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_8",
                        &BToKstarDilepton::integrated_j_8,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::J_9",
                        &BToKstarDilepton::integrated_j_9,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_3",
                        &BToKstarDilepton::integrated_s_3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_4",
                        &BToKstarDilepton::integrated_s_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_5",
                        &BToKstarDilepton::integrated_s_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_7",
                        &BToKstarDilepton::integrated_s_7,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_8",
                        &BToKstarDilepton::integrated_s_8,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::S_9",
                        &BToKstarDilepton::integrated_s_9,
                        std::make_tuple("q2_min", "q2_max")),

                // here we put observables S_X in the LHCb angular conveniton:
                make_observable("B->K^*ll::S_1s_LHCb",
                        &BToKstarDilepton::integrated_s_1s_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_1c_LHCb",
                        &BToKstarDilepton::integrated_s_1c_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_2s_LHCb",
                        &BToKstarDilepton::integrated_s_2s_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_2c_LHCb",
                        &BToKstarDilepton::integrated_s_2c_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_3_LHCb",
                        &BToKstarDilepton::integrated_s_3_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_4_LHCb",
                        &BToKstarDilepton::integrated_s_4_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_5_LHCb",
                        &BToKstarDilepton::integrated_s_5_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_6s_LHCb",
                        &BToKstarDilepton::integrated_s_6s_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_6c_LHCb",
                        &BToKstarDilepton::integrated_s_6c_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_7_LHCb",
                        &BToKstarDilepton::integrated_s_7_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_8_LHCb",
                        &BToKstarDilepton::integrated_s_8_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::S_9_LHCb",
                        &BToKstarDilepton::integrated_s_9_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::A_FB_LHCb",
                        &BToKstarDilepton::integrated_forward_backward_asymmetry_LHCb,
                        std::make_tuple("s_min", "s_max")),

                make_observable("B->K^*ll::A_9",
                        &BToKstarDilepton::integrated_a_9,
                        std::make_tuple("q2_min", "q2_max")),

#if 0
                make_observable("B->K^*ll::D_4",
                        &BToKstarDilepton::integrated_d_4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::D_5",
                        &BToKstarDilepton::integrated_d_5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("B->K^*ll::D_6s",
                        &BToKstarDilepton::integrated_d_6s,
                        std::make_tuple("q2_min", "q2_max")),
#endif

                make_observable_ratio("B->K^*ll::R_K^*",
                        R"(R_{K^*})",
                        &BToKstarDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_mu_min", "q2_mu_max"),
                        Options{ { "l", "mu" } },
                        &BToKstarDilepton::integrated_branching_ratio,
                        std::make_tuple("q2_e_min", "q2_e_max"),
                        Options{ { "l", "e" } }
                        ),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Lambda_b -> Lambda l^+ l^-
    // {{{
    ObservableGroup
    make_lambdab_to_lambda_ll_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $\Lambda_b \to \Lambda\ell^+\ell^-$ decays)",
            R"(The option "l" selects the charged lepton flavour.)",
            {
                // Lambda_b -> Lambda l^+ l^-, Large Recoil
                make_observable("Lambda_b->Lambdall::dBR/dq2@LargeRecoil", R"(d\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-)/dq^2)",
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^l(q2)@LargeRecoil",
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^h(q2)@LargeRecoil",
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^c(q2)@LargeRecoil",
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::F_0(q2)@LargeRecoil",
                        &LambdaBToLambdaDilepton<LargeRecoil>::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::BR@LargeRecoil", R"(\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^l@LargeRecoil", R"(A_\text{FB}^\ell(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^h@LargeRecoil", R"(A_\text{FB}^h(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^c@LargeRecoil", R"(A_\text{FB}^{h,\ell}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::F_0@LargeRecoil", R"(F_0(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LargeRecoil>::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                // Lambda_b -> Lambda l^+ l^-, Low Recoil
                make_observable("Lambda_b->Lambdall::dBR/dq2@LowRecoil", R"(d\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-)/dq^2)",
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^l(q2)@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_leptonic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^h(q2)@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_hadronic,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::A_FB^c(q2)@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_a_fb_combined,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::F_0(q2)@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::differential_fzero,
                        std::make_tuple("q2")),

                make_observable("Lambda_b->Lambdall::BR@LowRecoil", R"(\mathcal{B}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^l@LowRecoil", R"(A_\text{FB}^\ell(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_leptonic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^h@LowRecoil", R"(A_\text{FB}^h(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_hadronic,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::A_FB^c@LowRecoil", R"(A_\text{FB}^{h,\ell}(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_a_fb_combined,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::F_0@LowRecoil", R"(F_0(\Lambda_b\to\Lambda\ell^+\ell^-))",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_fzero,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1ss@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1cc@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_1c@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k1c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2ss@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2ss,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2cc@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2cc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_2c@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k2c,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_3sc@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_3s@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k3s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_4sc@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4sc,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::K_4s@LowRecoil",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_k4s,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_1@LowRecoil", R"(M_1)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m1,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_2@LowRecoil", R"(M_2)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m2,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_3@LowRecoil", R"(M_3)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m3,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_4@LowRecoil", R"(M_4)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m4,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_5@LowRecoil", R"(M_5)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m5,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_6@LowRecoil", R"(M_6)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m6,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_7@LowRecoil", R"(M_7)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m7,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_8@LowRecoil", R"(M_8)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m8,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_9@LowRecoil", R"(M_9)",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m9,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_10@LowRecoil", R"(M_{10})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m10,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_11@LowRecoil", R"(M_{11})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m11,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_12@LowRecoil", R"(M_{12})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m12,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_13@LowRecoil", R"(M_{13})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m13,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_14@LowRecoil", R"(M_{14})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m14,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_15@LowRecoil", R"(M_{15})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m15,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_16@LowRecoil", R"(M_{16})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m16,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_17@LowRecoil", R"(M_{17})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m17,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_18@LowRecoil", R"(M_{18})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m18,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_19@LowRecoil", R"(M_{19})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m19,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_20@LowRecoil", R"(M_{20})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m20,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_21@LowRecoil", R"(M_{21})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m21,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_22@LowRecoil", R"(M_{22})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m22,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_23@LowRecoil", R"(M_{23})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m23,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_24@LowRecoil", R"(M_{24})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m24,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_25@LowRecoil", R"(M_{25})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m25,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_26@LowRecoil", R"(M_{26})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m26,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_27@LowRecoil", R"(M_{27})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m27,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_28@LowRecoil", R"(M_{28})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m28,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_29@LowRecoil", R"(M_{29})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m29,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_30@LowRecoil", R"(M_{30})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m30,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_31@LowRecoil", R"(M_{31})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m31,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_32@LowRecoil", R"(M_{32})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m32,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_33@LowRecoil", R"(M_{33})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m33,
                        std::make_tuple("q2_min", "q2_max")),

                make_observable("Lambda_b->Lambdall::M_34@LowRecoil", R"(M_{34})",
                        &LambdaBToLambdaDilepton<LowRecoil>::integrated_m34,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // B -> X_s {gamma, l^+ l^-}
    // {{{
    ObservableGroup
    make_b_to_xs_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            R"(Observables in $B \to X_s \lbrace \gamma, \ell^+\ell^-\rbrace$ decays)",
            R"(The option "l" selects the charged lepton flavour. The option "q" selects the spectator quark flavour.)",
            {
                // B->X_s gamma
                make_observable("B->X_sgamma::BR@Minimal",
                        &BToXsGamma<Minimal>::integrated_branching_ratio),

                // B->X_s gamma, NLO implementation
                make_observable("B->X_sgamma::BR(E_min)@NLO",
                        &BToXsGamma<NLO>::integrated_branching_ratio,
                        std::make_tuple("E_min")),

                make_observable("B->X_sgamma::E_1(E_min)@NLO",
                        &BToXsGamma<NLO>::photon_energy_moment_1,
                        std::make_tuple("E_min")),

                make_observable("B->X_sgamma::E_2(E_min)@NLO",
                        &BToXsGamma<NLO>::photon_energy_moment_2,
                        std::make_tuple("E_min")),

                // B->X_s ll, HLMW2005
                make_observable("B->X_sll::dBR/dq2@HLMW2005",
                        &BToXsDilepton<HLMW2005>::differential_branching_ratio,
                        std::make_tuple("q2")),

                make_observable("B->X_sll::BR@HLMW2005",
                        &BToXsDilepton<HLMW2005>::integrated_branching_ratio,
                        std::make_tuple("q2_min", "q2_max")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    // Pseudo-observables for Non-local Matrix Elements
    // {{{
    ObservableGroup
    make_b_to_s_nonlocal_group()
    {
        auto imp = new Implementation<ObservableGroup>(
            "Pseudo-observables for Non-local Matrix Elements",
            "",
            {
                // B->K
                make_observable("B->Kgamma^*::Re{H_+}(q2)",
                        &NonlocalCorrelatorObservable<nc::BToK, nc::PToP>::re_H_plus,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B->Kgamma^*::Re{H_+}[s^1/s^0](q2)",
                        &NonlocalCorrelatorObservable<nc::BToK, nc::PToP>::re_normalized_moment_A,
                        std::make_tuple("q2")),

                // B->K^*
                make_observable("B->K^*gamma^*::Re{H_perp}(q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_H_perp,
                        std::make_tuple("q2")),

                make_observable("B->K^*gamma^*::Re{H_para}(q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_H_para,
                        std::make_tuple("q2")),

                make_observable("B->K^*gamma^*::Re{H_long}(q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_H_long,
                        std::make_tuple("q2")),

                // auxiliary variables, e.g. for determining the B-LCSR threshold parameters
                make_observable("B->K^*gamma^*::Re{H_1}[s^1/s^0](q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_normalized_moment_V1,
                        std::make_tuple("q2")),

                make_observable("B->K^*gamma^*::Re{H_2}[s^1/s^0](q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_normalized_moment_V2,
                        std::make_tuple("q2")),

                make_observable("B->K^*gamma^*::Re{H_23}[s^1/s^0](q2)",
                        &NonlocalCorrelatorObservable<nc::BToKstar, nc::PToV>::re_normalized_moment_V23,
                        std::make_tuple("q2")),
            }
        );

        return ObservableGroup(imp);
    }
    // }}}

    ObservableSection
    make_rare_b_decays_section()
    {
        auto imp = new Implementation<ObservableSection>(
            "Observables in rare (semi)leptonic and radiative $b$-hadron decays",
            "",
            {
                // B_q -> l^+ l^-
                make_b_to_ll_group(),

                // B_q -> P l^+ l^-
#if 0
                make_b_to_p_ll_group(),
#endif

                // B_q -> V gamma
                make_b_to_v_gamma_group(),

                // B_q -> V l^+ l^-
                make_b_to_v_ll_group(),

                // Lambda_b -> Lambda l^+ l^-
                make_lambdab_to_lambda_ll_group(),

                // B -> X_s {gamma, l^+ l^-}
                make_b_to_xs_group(),

                // Pseudo-observables for Non-local Matrix Elements
                make_b_to_s_nonlocal_group(),
            }
        );

        return ObservableSection(imp);
    }
}
