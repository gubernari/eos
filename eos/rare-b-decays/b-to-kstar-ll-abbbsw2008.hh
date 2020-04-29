/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_ABBBSW2008_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_LL_ABBBSW2008_HH 1

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-bfs2004.hh>

namespace eos
{
    template <>
    class BToKstarDileptonAmplitudes<tag::ABBBSW2008> :
        public BToKstarDileptonAmplitudes<tag::BFS2004>
    {
        public:
            BToKstarDileptonAmplitudes(const Parameters & p, const Options & o);
            ~BToKstarDileptonAmplitudes();

            virtual BToKstarDilepton::Amplitudes amplitudes(const double & q2) const;

            BToKstarDilepton::DipoleFormFactors dipole_form_factors(const double & q2, const WilsonCoefficients<BToS> & wc) const;
    };
}

#endif
