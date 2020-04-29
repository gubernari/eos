/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GP2004_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_GP2004_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::GP2004> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_c_MSbar;
            UsedParameter m_s;

            UsedParameter lambda_long;
            UsedParameter lambda_par;
            UsedParameter lambda_perp;

            UsedParameter sl_phase_long;
            UsedParameter sl_phase_par;
            UsedParameter sl_phase_perp;

            bool ccbar_resonance;

            bool use_nlo;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes();

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            inline complex<double> c7eff(const WilsonCoefficients<BToS> & wc, const double & q2) const;
            inline complex<double> c9eff(const WilsonCoefficients<BToS> & wc, const double & q2) const;
            inline double m_b_PS() const;
            inline double kappa() const;
            inline double norm(const double & q2) const;
    };
}

#endif
