/* vim: set sw=4 sts=4 tw=120 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Nico Gubernari
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

#ifndef MASTER_GUARD_EOS_FORM_FACTORS_BARYONIC_HQET_HH
#define MASTER_GUARD_EOS_FORM_FACTORS_BARYONIC_HQET_HH 1

#include <eos/form-factors/baryonic.hh>
#include <eos/utils/derivative.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/options.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/polylog.hh>
#include <eos/utils/power_of.hh>

#include <cmath>
#include <limits>

#include <iostream>

namespace eos
{
    using std::sqrt;

    /* HQET Form Factors, based on [BLPR2017] and [JS2018] */
    template <typename Process_, typename Transition_> class HQETFormFactors;

    class HQETVariables :
        public virtual ParameterUser
    {
      
        protected:
            std::shared_ptr<Model> _model;

            // spin avaraged baryon masses
            UsedParameter _m_Lb, _m_Lc;

            // parameter for modifying the z function
            UsedParameter _a;

            // option to determine the model for the leading-power IW function
            SwitchOption _opt_lp_model;
            std::function<double (const double &)> _zeta;

            // option to determine if we use z^3 terms in the leading-power IW function
            SwitchOption _opt_lp_zorder;
            double _enable_lp_z3;
            double _enable_lp_z4;
            double _enable_lp_z5;

            // option to determine if we use z^2 terms in the subleading-power IW function
            SwitchOption _opt_slp_zorder;
            double _enable_slp_z2;

            // option to determine if we use z^2 terms in the subsubleading-power IW function
            SwitchOption _opt_sslp_zorder;
            double _enable_sslp_z1;
            double _enable_sslp_z2;

            // option to determine if we use the SU3_F-symmetry limit for the subsubleading-power IW functions
            SwitchOption _opt_sslp_limit;

            // parameters for the leading Isgur-Wise function zeta
            UsedParameter _zetapone, _zetappone;




        public:
            HQETVariables(const Parameters & p, const Options & o) :
                _model(Model::make("SM", p, o)),
                _m_Lb(p["mass::Lambda_b"], *this),
                _m_Lc(p["mass::Lambda_c"], *this),
                _a(p["Lambda_b->Lambda_c::a@HQET"], *this),
                _opt_lp_model(o, "model-lp", { "power-series", "exponential" }, "power-series"),
                _opt_lp_zorder(o, "z-order-lp", { "2", "3", "4", "5" }, "3"),
                _enable_lp_z3(1.0 ? _opt_lp_zorder.value() >= "3" : 0.0),
                _enable_lp_z4(1.0 ? _opt_lp_zorder.value() >= "4" : 0.0),
                _enable_lp_z5(1.0 ? _opt_lp_zorder.value() >= "5" : 0.0),
                _opt_slp_zorder(o, "z-order-slp", { "1", "2" }, "2"),
                _enable_slp_z2(1.0 ? _opt_slp_zorder.value() >= "2" : 0.0),
                _opt_sslp_zorder(o, "z-order-sslp", { "0", "1", "2" }, "1"),
                _enable_sslp_z1(1.0 ? _opt_sslp_zorder.value() >= "1" : 0.0),
                _enable_sslp_z2(1.0 ? _opt_sslp_zorder.value() >= "2" : 0.0),
                _opt_sslp_limit(o, "SU3F-limit-sslp", { "0", "1" }, "0"),
                _zetapone(p["Lambda_b->Lambda_c::zeta'(1)@HQET"], *this),
                _zetappone(p["Lambda_b->Lambda_c::zeta''(1)@HQET"], *this)
            {
                using std::placeholders::_1;

                if (_opt_lp_model.value() == "exponential")
                {
                    _zeta = [=](const double & q2) -> double
                    {
                        return 0;
                    };
                }
                else
                {
                    _zeta = [=](const double & q2) -> double
                    {
                        return 0;
                    };
                }
            }



            ~HQETVariables() = default;


            virtual double _w(const double & q2) const = 0;

            double _q2(const double & w) const
            {
                return _m_Lb * _m_Lb + _m_Lc * _m_Lc -2.0 * _m_Lb * _m_Lc * w;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    results.add(Diagnostics::Entry{ _q2(1.0), "q2(w = 1.0)" });
                }

                return results;
            }
    };
}


#endif