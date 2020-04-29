/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-cffmpsv2015.hh>
#include <eos/rare-b-decays/nonlocal-correlator.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::CFFMPSV2015> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_c;
            UsedParameter m_s_MSbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            UsedParameter f_B;
            UsedParameter f_Kstar_par;
            UsedParameter f_Kstar_perp;
            UsedParameter lambda_B_p;
            UsedParameter a_1_par;
            UsedParameter a_2_par;
            UsedParameter a_1_perp;
            UsedParameter a_2_perp;

            UsedParameter uncertainty_xi_perp;
            UsedParameter uncertainty_xi_par;

            double e_q;

            char q;

            std::string coordinates;

            NonlocalCorrelatorPtr nonlocal_correlator;

            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &)> qcdf_dilepton_massless_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_charm_case;
            std::function<QCDFIntegrals<BToKstarDilepton> (const double &, const double &,
                    const double &, const double &, const double &, const double &,
                    const double &, const double &, const double &)> qcdf_dilepton_bottom_case;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes();

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            BToKstarDilepton::Amplitudes amp_semileptonic(const double & s, const WilsonCoefficients<BToS> & wc) const;
            WilsonCoefficients<BToS> wilson_coefficients() const;
            double m_b_PS() const;
            double mu_f() const;
            BToKstarDilepton::DipoleFormFactors dipole_form_factors(const double & q2, const WilsonCoefficients<BToS> & wc) const;
            double norm(const double & q2) const;
            double s_hat(const double & q2) const;
            double xi_perp(const double & q2) const;
            double xi_par(const double & q2) const;
    };
}

#endif
