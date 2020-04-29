/* vim: set sw=4 sts=4 et foldmethod=syntax : */

#ifndef MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BASE_HH
#define MASTER_GUARD_EOS_RARE_B_DECAYS_B_TO_KSTAR_GAMMA_BASE_HH 1

#include <eos/form-factors/mesonic.hh>
#include <eos/rare-b-decays/b-to-kstar-gamma.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options-impl.hh>

namespace eos
{
    class BToKstarGamma::AmplitudeGenerator :
        public ParameterUser
    {
        public:
            std::shared_ptr<Model> model;
            std::shared_ptr<FormFactors<PToV>> form_factors;

            UsedParameter hbar;

            UsedParameter mu;
            UsedParameter alpha_e;
            UsedParameter g_fermi;
            UsedParameter tau;

            UsedParameter m_B;
            UsedParameter m_Kstar;

            SwitchOption l;
            UsedParameter m_l;

            bool cp_conjugate;
            char q;
            double e_q;

            AmplitudeGenerator(const Parameters &, const Options &);

            virtual ~AmplitudeGenerator();
            virtual BToKstarGamma::Amplitudes amplitudes() const = 0;
    };

    template <typename Tag_> class BToKstarGammaAmplitudes;

    namespace tag
    {
        struct BCvDV2016;
        struct BFS2004;
    }
}

#endif
