/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Danny van Dyk
 * Copyright (c) 2011 Christian Wacker
 * Copyright (c) 2014 Frederik Beaujean
 * Copyright (c) 2014 Christoph Bobeth
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

#include <eos/rare-b-decays/b-to-kstar-ll-cffmpsv2015.hh>
#include <eos/rare-b-decays/charm-loops.hh>
#include <eos/rare-b-decays/hard-scattering.hh>
#include <eos/rare-b-decays/long-distance.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/destringify.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/memoise.hh>

#include <gsl/gsl_sf.h>

using namespace std;

namespace eos
{
    using namespace std::placeholders;
    using std::norm;

    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::BToKstarDileptonAmplitudes(const Parameters & p, const Options & o) :
        AmplitudeGenerator(p, o),
        hbar(p["hbar"], *this),
        m_b_MSbar(p["mass::b(MSbar)"], *this),
        m_c(p["mass::c"], *this),
        m_s_MSbar(p["mass::s(2GeV)"], *this),
        mu(p["mu"], *this),
        alpha_e(p["QED::alpha_e(m_b)"], *this),
        g_fermi(p["G_Fermi"], *this),
        tau(p["life_time::B_" + o.get("q", "d")], *this),
        f_B(p["decay-constant::B_" + o.get("q", "d")], *this),
        f_Kstar_par(p["B->K^*::f_Kstar_par"], *this),
        f_Kstar_perp(p["B->K^*::f_Kstar_perp@2GeV"], *this),
        lambda_B_p(p["lambda_B_p"], *this),
        a_1_par(p["B->K^*::a_1_par"], *this),
        a_2_par(p["B->K^*::a_2_par"], *this),
        a_1_perp(p["B->K^*::a_1_perp"], *this),
        a_2_perp(p["B->K^*::a_2_perp"], *this),
        uncertainty_xi_perp(p["formfactors::xi_perp_uncertainty"], *this),
        uncertainty_xi_par(p["formfactors::xi_par_uncertainty"], *this),
        e_q(-1.0/3.0)
    {
        std::string spectator_quark = o.get("q", "d");
        if ("d" == spectator_quark)
        {
            e_q = -1.0 / 3.0;
        }
        else if ("u" == spectator_quark)
        {
            e_q = 2.0 / 3.0;
        }
        else
        {
            throw InvalidOptionValueError("q", spectator_quark, "u, d");
        }

        // Select the appropriate parametrization of the non-factorizing
        // corrections at next-to-leading power
        std::string correlator = o.get("correlator", "CFFMPSV2015-cartesian");
        nonlocal_correlator = NonlocalCorrelator::make(correlator, p, o);
        if (! nonlocal_correlator)
        {
            throw InvalidOptionValueError("correlator", correlator, "CFFMPSV2015-cartesian, CFFMPSV2015-polar, BCvD2016-model1-cartesian, BCvD2016-model1-polar");
        }

        // Select the appropriate calculator for the QCDF integrals
        std::string qcdf_integrals(o.get("qcdf-integrals", "mixed"));
        if ("mixed" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Mixed>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else if ("numerical" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Numerical>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else if ("analytical" == qcdf_integrals)
        {
            qcdf_dilepton_massless_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_massless_case,
                        _1, _2, _3, _4, _5, _6, _7, _8);
            qcdf_dilepton_charm_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_charm_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
            qcdf_dilepton_bottom_case = std::bind(&QCDFIntegralCalculator<BToKstarDilepton, tag::Analytical>::dilepton_bottom_case,
                        _1, _2, _3, _4, _5, _6, _7, _8, _9);
        }
        else
        {
            throw InvalidOptionValueError("qcdf-integrals", qcdf_integrals, "mixed, numerical, analytical");
        }
    }

    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::~BToKstarDileptonAmplitudes()
    {
    }

    WilsonCoefficients<BToS>
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::wilson_coefficients() const
    {
        return model->wilson_coefficients_b_to_s(lepton_flavour, cp_conjugate);
    }

    // as in ABBBSW2008
    BToKstarDilepton::DipoleFormFactors
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::dipole_form_factors(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        // charges of down- and up-type quarks
        static const double
            e_d = -1.0 / 3.0,
            e_u = +2.0 / 3.0;

        // spectator contributions
        const double delta_qu = (q == 'u' ? 1.0 : 0.0);

        // kinematics
        const double
            m_c_pole = model->m_c_pole(),
            m_b_PS = this->m_b_PS(),
            energy = this->energy(s);

        // couplings
        const double
            alpha_s_mu = model->alpha_s(mu()), // alpha_s at the hard scale
            a_mu = alpha_s_mu * QCD::casimir_f / 4.0 / M_PI,
            alpha_s_mu_f = model->alpha_s(std::sqrt(mu() * 0.5)), // alpha_s at the factorization scale
            a_mu_f = alpha_s_mu_f * QCD::casimir_f / 4.0 / M_PI;

        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);

        QCDFIntegrals<BToKstarDilepton>
            qcdf_0 = this->qcdf_dilepton_massless_case(s, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par),
            qcdf_c = this->qcdf_dilepton_charm_case(s, m_c_pole, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par),
            qcdf_b = this->qcdf_dilepton_bottom_case(s, m_b_PS, m_B, m_Kstar, mu, a_1_perp, a_2_perp, a_1_par, a_2_par);

        // inverse of the "negative" moment of the B meson LCDA
        // cf. [BFS2001], Eq. (54), p. 15
        const double
            omega_0 = lambda_B_p,
            lambda_B_p_inv = 1.0 / lambda_B_p;

        const complex<double>
            lambda_B_m_inv = complex<double>(-gsl_sf_expint_Ei(s / m_B / omega_0), M_PI) * (std::exp(-s / m_B / omega_0) / omega_0);

        /* Effective wilson coefficients */
        // cf. [BFS2001], below Eq. (26), p. 8
        const complex<double> c8eff = wc.c8() + wc.c3() - 1.0/6.0 * wc.c4() + 20.0 * wc.c5() - 10.0/3.0 * wc.c6();

        /* perpendicular, top sector */
        // cf. [BFS2001], Eqs. (34), (37), p. 9
        const complex<double>
            C1nf_top_perp = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (s / (2.0 * m_b_PS * m_B)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                    + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS))),

        /* perpendicular, up sector */
        // cf. [BFS2001], Eqs. (34), (37), p. 9
        // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            C1nf_up_perp = (-1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (s / (2.0 * m_b_PS * m_B)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS)))),

        /* parallel, top sector */
        // cf. [BFS2001], Eqs. (38), p. 9
            C1nf_top_par = (+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole)
                + c8eff * CharmLoops::F87_massless(mu, s, m_b_PS)
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole)
                    + wc.c2() * memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole)
                    + c8eff * CharmLoops::F89_massless(s, m_b_PS))),

        /* parallel, up sector */
        // cf. [BFS2004], last paragraph in Sec A.1, p. 24
        // [BFS2004], [S2004] have a different sign convention for F{12}{79}_massless than we!
            C1nf_up_par = (+1.0 / QCD::casimir_f) * (
                (wc.c2() - wc.c1() / 6.0) * (memoise(CharmLoops::F27_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F27_massless(mu, s, m_b_PS))
                + (m_B / (2.0 * m_b_PS)) * (
                    wc.c1() * (memoise(CharmLoops::F19_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F19_massless(mu, s, m_b_PS))
                    + wc.c2() * (memoise(CharmLoops::F29_massive, mu(), s, m_b_PS, m_c_pole) - CharmLoops::F29_massless(mu, s, m_b_PS))));

        // compute the factorizing contributions
        // in ABBBSW2008: C0 is included in naively factorizing part and C1f = 0
        const complex<double>
            C_perp = a_mu * (C1nf_top_perp + lambda_hat_u * C1nf_up_perp),
            C_par  = a_mu * (C1nf_top_par  + lambda_hat_u * C1nf_up_par);


        /* perpendicular, top sector */
        // cf. [BFS2001], Eq. (23), p. 7
        const complex<double>
            T1nf_top_perp_p = (-4.0 * e_d * c8eff * qcdf_0.j0bar_perp
                + m_B / (2.0 * m_b_PS) * (
                    e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0/3.0 * wc.c6() - (4.0 * m_b_PS / m_B) * (wc.c3() - wc.c4()/6.0 + 4.0 * wc.c5() - 2.0/3.0 * wc.c6())) * qcdf_b.jtilde1_perp
                    + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() - 8.0/3.0 * wc.c6()) * qcdf_0.jtilde1_perp)) * lambda_B_p_inv,
        // T1nf_top_perp_m = 0, cf. [BFS2001], Eq. (17), p. 6

        /* perpendicular, up sector */
        // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
        // cf. [BFS2004], Eq. (50), p. 25
            T1nf_up_perp_p = +e_u * m_B / (2.0 * m_b_PS) * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde1_perp - qcdf_0.jtilde1_perp) * lambda_B_p_inv,

        /* parallel, top sector */
        // cf. [BFS2004], Eqs. (46)-(47), p. 25 without the \omega term.
            T0_top_par_m = -e_q * 4.0 * m_B / m_b_PS * (wc.c3() + 4.0/3.0 * wc.c4() + 16.0 * wc.c5() + 64.0/3.0 * wc.c6()) * lambda_B_m_inv,
        // cf. [BFS2001], Eq. (25), p. 7
            T1nf_top_par_p = m_B / m_b_PS * (
                e_u * (-wc.c1() / 6.0 + wc.c2() + 6.0 * wc.c6()) * qcdf_c.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() + 10.0 / 3.0 * wc.c6()) * qcdf_b.jtilde2_parallel
                + e_d * (wc.c3() - wc.c4() / 6.0 + 16.0 * wc.c5() -  8.0 / 3.0 * wc.c6()) * qcdf_0.jtilde2_parallel) * lambda_B_p_inv,
        // cf. [BFS2001], Eq. (26), pp. 7-8
            T1nf_top_par_m = e_q * (8.0 * c8eff * qcdf_0.j0_parallel
                + 6.0 * m_B / m_b_PS * (
                    (-wc.c1() / 6.0 + wc.c2() + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j4_parallel
                    + (wc.c3() + 5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j4_parallel
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j4_parallel
                    -8.0 / 27.0 * (-7.5 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))) * lambda_B_m_inv;

        /* parallel, up sector */
        // all T1f_up vanish, cf. [BFS2004], sentence below Eq. (49), p. 25
        // cf. [BFS2004], Eqs. (46),(48), p. 25 without the \omega term
        const complex<double>
            T0_up_par_m = +e_q * 4.0 * m_B / m_b_PS * (3.0 * delta_qu * wc.c2()) * lambda_B_m_inv,
        // cf. [BFS2004], Eq. (50), p. 25
            T1nf_up_par_p = +e_u * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.jtilde2_parallel - qcdf_0.jtilde2_parallel) * lambda_B_p_inv,
        // cf. [BFS2004], Eq. (50), p. 25 without the \omega term
            T1nf_up_par_m = +e_q * 6.0 * m_B / m_b_PS * (-wc.c1() / 6.0 + wc.c2()) * (qcdf_c.j4_parallel - qcdf_0.j4_parallel) * lambda_B_m_inv;


        // Compute the nonfactorizing contributions
        const complex<double>
            T_perp = a_mu_f * (T1nf_top_perp_p + lambda_hat_u * T1nf_up_perp_p),
            T_par  = a_mu_f * (T1nf_top_par_p  + lambda_hat_u * T1nf_up_par_p)
                     + T0_top_par_m + lambda_hat_u * T0_up_par_m
                     + a_mu_f * (T1nf_top_par_m + lambda_hat_u * T1nf_up_par_m);

        // Compute the numerically leading power-suppressed weak annihilation contributions to order alpha_s^0
        // cf. [BFS2004], Eq. (51)
        const complex<double>
            Delta_T_ann_top_perp = e_q * M_PI * M_PI * f_B / 3.0 / m_b_PS / m_B * (
                -4.0 * f_Kstar_perp * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 3.0 * wc.c5() + 4.0 * wc.c6())) * qcdf_0.j0_perp
                + 2.0 * f_Kstar_par * (wc.c3() + 4.0 / 3.0 * (wc.c4() + 12.0 * wc.c5() + 16.0 * wc.c6())) *
                    (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p)),
            Delta_T_ann_up_perp = -e_q * 2.0 * M_PI * M_PI * f_B * f_Kstar_par / 3.0 / m_b_PS / m_B
                * (m_Kstar / (1.0 - s / (m_B * m_B)) / lambda_B_p) * 3.0 * delta_qu * wc.c2(),
        // Compute the numerically leading power-suppressed hard spectator interaction contributions to order alpha_s^1
        // cf. [BFS2004], Eqs. (52), (53)
            Delta_T_hsa_top_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                12.0 * c8eff * (m_b_PS / m_B) * f_Kstar_perp() * 1.0 / 3.0 * (qcdf_0.j0_perp + qcdf_0.j7_perp)
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (
                      (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j5_perp
                    + (wc.c3() +  5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j5_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j5_perp
                    - (8.0 / 27.0) * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()) * qcdf_0.j0_perp)
                - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (
                      (wc.c2() - wc.c1() / 6.0 + wc.c4() + 10.0 * wc.c6()) * qcdf_c.j6_perp
                    + (wc.c3() +  5.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 22.0 / 3.0 * wc.c6()) * qcdf_b.j6_perp
                    + (wc.c3() + 17.0 / 6.0 * wc.c4() + 16.0 * wc.c5() + 82.0 / 3.0 * wc.c6()) * qcdf_0.j6_perp
                    - 8.0 / 27.0 * (-15.0 / 2.0 * wc.c4() + 12.0 * wc.c5() - 32.0 * wc.c6()))),
            Delta_T_hsa_up_perp = e_q * a_mu_f * (M_PI * M_PI * f_B / (3.0 * m_b_PS * m_B)) * (
                + 8.0 * f_Kstar_perp * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0) * (qcdf_c.j5_perp - qcdf_0.j5_perp)
                - (4.0 * m_Kstar * f_Kstar_par / (1.0 - s / (m_B * m_B)) / lambda_B_p) * (3.0 / 4.0) * (wc.c2() - wc.c1() / 6.0)
                    * (qcdf_c.j6_perp - qcdf_0.j6_perp));

        // Compute the sum of the numerically leading power-suppressed contributions
        const complex<double>
            Delta_T_top_perp = Delta_T_ann_top_perp + Delta_T_hsa_top_perp,
            Delta_T_up_perp  = Delta_T_ann_up_perp + Delta_T_hsa_up_perp,
            Delta_T_perp     = Delta_T_top_perp + lambda_hat_u * Delta_T_up_perp;

        // cf. [BFS2001], Eq. (15), and [BHP2008], Eq. (C.4)
        BToKstarDilepton::DipoleFormFactors result;
        result.calT_perp_left  = xi_perp(s) * C_perp + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_perp) / m_B * T_perp + Delta_T_perp;
        result.calT_perp_right = result.calT_perp_left;
        result.calT_parallel   = xi_par(s) * C_par   + power_of<2>(M_PI) / 3.0 * (f_B * f_Kstar_par * m_Kstar) / (m_B * energy) * T_par;

        return result;
    }

    /* Form factors */
    //  cf. [BHP2008], Eq. (E.4), p. 23
    double
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::xi_perp(const double & s) const
    {
        const double factor = m_B() / (m_B() + m_Kstar());
        double result = uncertainty_xi_perp * factor * form_factors->v(s);

        return result;
    }

    double
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::xi_par(const double & s) const
    {
        const double factor1 = (m_B() + m_Kstar()) / (2.0 * energy(s));
        const double factor2 = (1.0 - m_Kstar() / m_B());
        double result = uncertainty_xi_par * (factor1 * form_factors->a_1(s) - factor2 * form_factors->a_2(s));

        return result;
    }

    double BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::norm(const double & s) const
    {
        double lambda_t2 = std::norm(model->ckm_tb() * conj(model->ckm_ts()));

        return g_fermi() * alpha_e() * std::sqrt(
                  1.0 / 3.0 / 1024 / power_of<5>(M_PI) / m_B()
                  * lambda_t2 * s_hat(s) * std::sqrt(lambda(s)) * beta_l(s)
               ); // cf. [BHP2008], Eq. (C.6), p. 21
    }

    double
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::s_hat(const double & s) const
    {
        return s / m_B() / m_B();
    }

    double
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::mu_f() const
    {
        return 1.5;
    }

    double
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::m_b_PS() const
    {
        // Actually use the PS mass at mu_f = 1.5 GeV
        return model->m_b_ps(mu_f());
    }

    /* Amplitudes */
#if 0
    /* calculate contribution to transversity amplitudes
     * -> from semi-leptonic operators O_9,10,S,T + primed and
     *    el-mg. dipole operator O_7, which factorize naively
     * -> not included 4-quark matrix elements to C_7,9 -> C_7,9^eff,
     *    which are in hadronic contributions to transversity amplitudes
     * -> use full QCD form factors (FF)
     *    TODO: implement variants for large and low recoil FF-relations
     *          a la BFS2004 for parts ~ C_7,7'
     * -> buffer intermediate results for potential later use
     */
    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::amp_semileptonic(const double & s, const WilsonCoefficients<BToS> & wc) const
    {
        BToKstarDilepton::Amplitudes amps;

        const double
            shat = s_hat(s),
            sqrt_s = std::sqrt(s),
            mKhat2 = power_of<2>(m_Kstar() / m_B()),
            m_K2 = power_of<2>(m_Kstar()),
            m_B2 = power_of<2>(m_B()),
            m_sum = m_B() + m_Kstar(),
            m_diff = m_B() - m_Kstar(),
            m2_diff = m_B2 - m_K2,
            norm_s = this->norm(s),
            lam_s = lambda(s),
            sqrt_lam = std::sqrt(lam_s);

        const double
            ff_V  = form_factors->v(s),
            ff_A0 = form_factors->a_0(s),
            ff_A1 = form_factors->a_1(s),
            ff_A2 = form_factors->a_2(s),
            ff_T1 = form_factors->t_1(s),
            ff_T2 = form_factors->t_2(s),
            ff_T3 = form_factors->t_3(s);

#if 0
        const double
            m_b_PS = this->m_b_PS();
#endif

        const complex<double>
            c910_mi_r = (wc.c9() - wc.c9prime()) + (wc.c10() - wc.c10prime()),
            c910_mi_l = (wc.c9() - wc.c9prime()) - (wc.c10() - wc.c10prime()),
            c910_pl_r = (wc.c9() + wc.c9prime()) + (wc.c10() + wc.c10prime()),
            c910_pl_l = (wc.c9() + wc.c9prime()) - (wc.c10() + wc.c10prime()),
            // here b-quark mass in O_7 is MS-bar contrary to PS-scheme in [BFS2001]
            c7_mi = 2.0 * m_b_MSbar * (wc.c7() - wc.c7prime()),
            c7_pl = 2.0 * m_b_MSbar * (wc.c7() + wc.c7prime());

        // longitudinal, perpendicular, parallel amplitudes
        const double    // prefactors
            pre_long = -norm_s / (2.0 * m_Kstar() * sqrt_s),
            pre_perp = +std::sqrt(2.0) * norm_s * sqrt_lam,
            pre_para = -std::sqrt(2.0) * norm_s * m2_diff;

        const double
            kin_long_1 = (m2_diff - s) * m_sum * ff_A1 - lam_s / m_sum * ff_A2,
            kin_long_2 = (m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam_s / m2_diff * ff_T3;

        amps.a_long_right = pre_long * (c910_mi_r * kin_long_1 + c7_mi * kin_long_2);
        amps.a_long_left  = pre_long * (c910_mi_l * kin_long_1 + c7_mi * kin_long_2);

        amps.a_perp_right = pre_perp * (c910_pl_r * ff_V / m_sum + c7_pl / s * ff_T1);
        amps.a_perp_left  = pre_perp * (c910_pl_l * ff_V / m_sum + c7_pl / s * ff_T1);

        amps.a_para_right = pre_para * (c910_mi_r * ff_A1 / m_diff + c7_mi / s * ff_T2);
        amps.a_para_left  = pre_para * (c910_mi_l * ff_A1 / m_diff + c7_mi / s * ff_T2);

        // timelike amplitude
        amps.a_time = norm_s * sqrt_lam / sqrt_s
            * (2.0 * (wc.c10() - wc.c10prime())
               + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime())
              ) * ff_A0;

        // scalar amplitude
        amps.a_scal = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * ff_A0;

        // tensor amplitudes
        const double
            kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam_s / m2_diff * ff_T3),
            kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
            kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

        // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
        static const double sign = -1;

        amps.a_para_perp = kin_tensor_1 * wc.cT();
        amps.a_time_long = kin_tensor_1 * sign * wc.cT5();

        amps.a_time_perp = kin_tensor_2 * wc.cT();
        amps.a_long_perp = kin_tensor_2 * sign * wc.cT5();

        amps.a_time_para = kin_tensor_3 * sign * wc.cT5();
        amps.a_long_para = kin_tensor_3 * wc.cT();

        // buffer intermediate results for potential later use
        cbuf.shat = shat; cbuf.sqrt_s = sqrt_s; cbuf.mKhat2 = mKhat2;
        cbuf.m_K2 = m_K2; cbuf.m_B2 = m_B2;
        cbuf.m_sum = m_sum; cbuf.m_diff = m_diff; cbuf.m2_diff = m2_diff;
        cbuf.norm_s = norm_s; cbuf.lam_s = lam_s; cbuf.sqrt_lam = sqrt_lam;

        cbuf.ff_V = ff_V;
        cbuf.ff_A0 = ff_A0; cbuf.ff_A1 = ff_A1; cbuf.ff_A2 = ff_A2;
        cbuf.ff_T1 = ff_T1; cbuf.ff_T2 = ff_T2; cbuf.ff_T3 = ff_T3;

        cbuf.pre_long = pre_long; cbuf.pre_perp = pre_perp; cbuf.pre_par = pre_para;
        cbuf.kin_long_1 = kin_long_1; cbuf.kin_long_2 = kin_long_2;

        cbuf.c910_mi_r = c910_mi_r; cbuf.c910_mi_l = c910_mi_l;
        cbuf.c910_pl_r = c910_pl_r; cbuf.c910_pl_l = c910_pl_l;
        cbuf.c7_mi = c7_mi; cbuf.c7_pl = c7_pl;

        return amps;
    }
#endif

    // cf. [ABBBSW2008] and [CFFMPSV2015], use [BHvD2012] for tensor amplitudes
    // use full QCD form factors in leading QCDF (naively factorizing) amplitudes
    // use soft form factors in non-factorizable contributions (~ alpha_s)
    BToKstarDilepton::Amplitudes
    BToKstarDileptonAmplitudes<tag::CFFMPSV2015>::amplitudes(const double & s) const
    {
        BToKstarDilepton::Amplitudes result;

        WilsonCoefficients<BToS> wc = wilson_coefficients();

        const double
            sqrt_s = std::sqrt(s),
            m_K2 = power_of<2>(m_Kstar()),
            m_B2 = power_of<2>(m_B()),
            m_sum = m_B() + m_Kstar(),
            m_diff = m_B() - m_Kstar(),
            m2_diff = m_B2 - m_K2,
            norm_s = this->norm(s),
            lam_s = lambda(s),
            sqrt_lam = std::sqrt(lam_s);

        // full form factors:
        const double
            ff_V   = form_factors->v(s),
            ff_A0  = form_factors->a_0(s),
            ff_A1  = form_factors->a_1(s),
            ff_A2  = form_factors->a_2(s),
            ff_T1  = form_factors->t_1(s),
            ff_T2  = form_factors->t_2(s),
            ff_T3  = form_factors->t_3(s);

        // quark masses:
        const double
            m_c_pole = model->m_c_pole(),
            m_b_PS = this->m_b_PS();

        // charge conjugate
        complex<double> lambda_hat_u = (model->ckm_ub() * conj(model->ckm_us())) / (model->ckm_tb() * conj(model->ckm_ts()));
        if (cp_conjugate)
            lambda_hat_u = std::conj(lambda_hat_u);

        /* Y(s) for the up and the top sector for effective Wilson coefficients */
        // cf. [BFS2001s], Eq. (10), p. 4
        // This is the corretion to SM WCs.
        const complex<double>
            Y_top_c = 4.0 / 3.0 * wc.c1() + wc.c2() + 6.0 * wc.c3() + 60.0 * wc.c5(),
            Y_top_b = -0.5 * (7.0 * wc.c3() + 4.0 / 3.0 * wc.c4() + 76.0 * wc.c5() + 64.0 / 3.0 * wc.c6()),
            Y_top_0 = -0.5 * (wc.c3() + 4.0 / 3.0 * wc.c4() + 16.0 * wc.c5() + 64 / 3.0 * wc.c6()),
            Y_top_  = 2.0 / 9.0 * (6.0 * wc.c3() + 32.0 * wc.c5() + 32.0 / 3.0 * wc.c6());

        // Use b pole mass according to [BFS2001], Sec. 3.1, paragraph Quark Masses,
        // then replace b pole mass by the PS mass.
        complex<double> Y_top = Y_top_c * CharmLoops::h(mu, s, m_c_pole)
             + Y_top_b * CharmLoops::h(mu, s, m_b_PS)
             + Y_top_0 * CharmLoops::h(mu, s)
             + Y_top_;
        // cf. [BFS2004], Eq. (43), p. 24
        complex<double> Y_up = (4.0 / 3.0 * wc.c1() + wc.c2()) * (CharmLoops::h(mu, s, m_c_pole) - CharmLoops::h(mu, s));

        // calculating effective wilson coeficients
        const complex<double>
            // cf. [BFS2001], below Eq. (9), p. 4
            c7eff = - 1.0/3.0 * wc.c3() - 4.0/9.0 * wc.c4() - 20.0/3.0 * wc.c5() - 80.0/9.0 * wc.c6(), // modification to Delta C7, aka Delta C7
            c9eff = Y_top + lambda_hat_u * Y_up,  // modification to C9
            c910_mi_r = (wc.c9() + c9eff - wc.c9prime()) + (wc.c10() - wc.c10prime()),
            c910_mi_l = (wc.c9() + c9eff - wc.c9prime()) - (wc.c10() - wc.c10prime()),
            c910_pl_r = (wc.c9() + c9eff + wc.c9prime()) + (wc.c10() + wc.c10prime()),
            c910_pl_l = (wc.c9() + c9eff + wc.c9prime()) - (wc.c10() + wc.c10prime()),
            // here b-quark mass in O_7 is MS-bar contrary to PS-scheme in [BFS2001]
            c7_mi = 2.0 * m_b_MSbar * (wc.c7() + c7eff - wc.c7prime()),
            c7_pl = 2.0 * m_b_MSbar * (wc.c7() + c7eff + wc.c7prime());

        // cf. [CFFMPSV2015], eq.(2.7)
        // Note that [KMPW2010] use interchanged definition of color-octet and singlet charged-current operators O_1,2 (see App. A)
        // w.r.t. standard literature => Wilson coefficient C_1[KMPW2010] -> C_2[standard]
        // but this is irrelevant, because CFFMPSV2015 have multiplied g_M1,2,3 by the same factor and the Wilson coefficient cancels
        complex<double> h_0 = nonlocal_correlator->h_long(s);
        complex<double> h_p = nonlocal_correlator->h_plus(s);
        complex<double> h_m = nonlocal_correlator->h_minus(s);
        const complex<double> g_M_fac = m_B2 * m_B() * power_of<2>(4 * M_PI) / s;
        const complex<double> g_M1 = -1. * g_M_fac * m_sum/sqrt_lam * (h_m - h_p);
        const complex<double> g_M2 = -1. * g_M_fac / m_sum * (h_m + h_p);
        const complex<double> g_M3 =       g_M_fac * m_sum / lam_s
                                   * ( 4. * m_Kstar * sqrt_s * h_0 - (m2_diff - s) * (h_m + h_p));

        // longitudinal amplitude
        const double pre_long = -norm_s / (2.0 * m_Kstar() * sqrt_s);
        const complex<double>
          a = (m2_diff - s) * m_sum * ff_A1 - lam_s / m_sum * ff_A2,
          b = (m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam_s / m2_diff * ff_T3,
          c = (m2_diff - s) * m_sum * g_M2 - lam_s / m_sum * g_M3;

        result.a_long_right = pre_long * (c910_mi_r * a + c + c7_mi * b);
        result.a_long_left  = pre_long * (c910_mi_l * a + c + c7_mi * b);

        // perpendicular amplitude
        const double pre_perp = +std::sqrt(2.0) * norm_s * sqrt_lam;

        result.a_perp_right = pre_perp * ( (c910_pl_r * ff_V + g_M1) / m_sum + c7_pl / s * ff_T1);
        result.a_perp_left  = pre_perp * ( (c910_pl_l * ff_V + g_M1) / m_sum + c7_pl / s * ff_T1);

        // parallel amplitude
        const double pre_para = -std::sqrt(2.0) * norm_s * m2_diff;

        result.a_para_right = pre_para * ( (c910_mi_r * ff_A1 + g_M2) / m_diff + c7_mi / s * ff_T2);
        result.a_para_left  = pre_para * ( (c910_mi_l * ff_A1 + g_M2) / m_diff + c7_mi / s * ff_T2);

        // timelike amplitude
        result.a_time = norm_s * sqrt_lam / sqrt_s
            * (2.0 * (wc.c10() - wc.c10prime())
               + s / m_l / (m_b_MSbar + m_s_MSbar) * (wc.cP() - wc.cPprime())
              ) * ff_A0;

        // scalar amplitude
        result.a_scal = -2.0 * norm_s * sqrt_lam * (wc.cS() - wc.cSprime()) / (m_b_MSbar + m_s_MSbar) * ff_A0;

        // tensor amplitudes
        const double
            kin_tensor_1 = norm_s / m_Kstar() * ((m_B2 + 3.0 * m_K2 - s) * ff_T2 - lam_s / m2_diff * ff_T3),
            kin_tensor_2 = 2.0 * norm_s * sqrt_lam / sqrt_s * ff_T1,
            kin_tensor_3 = 2.0 * norm_s * m2_diff / sqrt_s * ff_T2;

        // correct the sign of C_T5 from [BHvD2012v4] because of inconsistent use of gamma5 <-> Levi-Civita
        static const double sign = -1;

        result.a_para_perp = kin_tensor_1 * wc.cT();
        result.a_time_long = kin_tensor_1 * sign * wc.cT5();

        result.a_time_perp = kin_tensor_2 * wc.cT();
        result.a_long_perp = kin_tensor_2 * sign * wc.cT5();

        result.a_time_para = kin_tensor_3 * sign * wc.cT5();
        result.a_long_para = kin_tensor_3 * wc.cT();

        //
        // Beyond Naive factorization part - from QCDF
        //

        auto dff = dipole_form_factors(s, wc);

        // these kinematical factors reduce for mKstar = 0 to [ABBBSW2008] eq. (3.46)
        // using approximation mKstar = 0 as in [ABBBSW2008] eq. (3.46)
        // [christoph] !!! setting mKstar = 0 causes quite large shifts in Br (+17%), F_L (-17%) (integrated from 1.1 to 6 GeV^2, KMPW-form factors)
        // did not test other observables

        // the Beyond naive factorization part is nadled by the CFFMPSV2015 aproach

        const double
            pre_p = std::sqrt(2.0) * norm_s * 2.0 * m_b_MSbar / s * (m_B2 - s),

            pre_0 = norm_s / m_Kstar() / m_B2 / sqrt_s * m_b_MSbar * power_of<2>(m_B2 - s);

        // here we assume that there is no nonfactorizable correction (hence 1. * )
        result.a_long_right += 1. * pre_0 * dff.calT_parallel;
        result.a_long_left  += 1. * pre_0 * dff.calT_parallel;

        result.a_perp_right += 1. * pre_p * dff.calT_perp_right;
        result.a_perp_left  += 1. * pre_p * dff.calT_perp_right;

        result.a_para_right -= 1. * pre_p * dff.calT_perp_left;
        result.a_para_left  -= 1. * pre_p * dff.calT_perp_left;

        return result;
    }
}
