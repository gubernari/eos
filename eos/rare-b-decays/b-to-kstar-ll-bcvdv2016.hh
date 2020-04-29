/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_CFFMPSV2015_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-cffmpsv2015.hh>
#include <eos/rare-b-decays/nonlocal-correlator.hh>
#include <eos/rare-b-decays/qcdf-integrals.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::BCvDV2016> :
        public BToKstarDilepton::AmplitudeGenerator
    {
        public:
            UsedParameter hbar;

            UsedParameter m_b_MSbar;
            UsedParameter m_s_MSbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            SwitchOption q;

            SwitchOption correlator;
            NonlocalCorrelatorPtr<nc::PToV> nonlocal_correlator;

            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes() = default;

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            WilsonCoefficients<BToS> wilson_coefficients() const;
    };
}

#endif
