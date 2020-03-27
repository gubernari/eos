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
            UsedParameter _zetapone, _zetappone, _zetapppone, _zetappppone, _zetapppppone;




        public:
            HQETVariables(const Parameters & p, const Options & o) :
                _model(Model::make("SM", p, o)),
                _m_Lb(p["mass::Lambda_b"], *this),
                _m_Lc(p["mass::Lambda_c"], *this),
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
                _zetappone(p["Lambda_b->Lambda_c::zeta''(1)@HQET"], *this),
                _zetapppone(p["Lambda_b->Lambda_c::zeta'''(1)@HQET"], *this),
                _zetappppone(p["Lambda_b->Lambda_c::zeta''''(1)@HQET"], *this),
                _zetapppppone(p["Lambda_b->Lambda_c::zeta'''''(1)@HQET"], *this)
            {}

            ~HQETVariables() = default;


            double _w(const double & q2) const
            {
                const double m_Lb = this->_m_Lb(), m_Lb2 = power_of<2>(m_Lb);
                const double m_Lc = this->_m_Lc(), m_Lc2 = power_of<2>(m_Lc);

                return (m_Lb2 + m_Lc2 - q2) / (2.0 * m_Lb * m_Lc);
            }

            double _q2(const double & w) const
            {
                return _m_Lb * _m_Lb + _m_Lc * _m_Lc -2.0 * _m_Lb * _m_Lc * w;
            }

            double _z(const double & w) const
            {
                return 0.0;
            }

            const double _sp() const
            {
                  return pow(this->_m_Lb() + this->_m_Lb(), 2);
            }

            const double _s0() const
            {
                  return pow(this->_m_Lb() - this->_m_Lb(), 2);
            }

            double _zeta_power_series(const double & q2) const
            {
                const double m_Lb = this->_m_Lb();
                const double m_Lc = this->_m_Lc();

                const double sp = this->_sp();
                const double s0 = this->_s0();

                // expansion in z around z_0
                const double  z_0 = _z(_q2(1.0));
                const double  z   = (_z(q2) - z_0);
                const double z2   =  z *  z;
                const double z3   = z2 *  z * _enable_lp_z3;
                const double z4   = z2 * z2 * _enable_lp_z4;
                const double z5   = z3 * z2 * _enable_lp_z5;


                const double wm11z1 = -(pow(m_Lb,-1)*pow(m_Lc,-1)*pow(-s0 + sp,-0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)*
                                      pow(-s0 + pow(m_Lb - m_Lc,2),2)*
                                      pow(s0 - 2*sp + pow(m_Lb - m_Lc,2) +
                                      2*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5),-1))/2.0;

                const double wm11z2 = (pow(m_Lb,-1)*pow(m_Lc,-1)*pow(s0 - sp,-1)*
                                      (14*s0*sp + 12*m_Lc*pow(m_Lb,3) - 3*pow(m_Lb,4) -
                                      3*pow(m_Lc,4) - pow(s0,2) - 16*pow(sp,2) +
                                      6*s0*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      16*sp*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      4*m_Lb*m_Lc*(6*s0 - 9*sp + 3*pow(m_Lc,2) -
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      2*pow(m_Lc,2)*(-6*s0 + 9*sp +
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      2*pow(m_Lb,2)*(-6*s0 + 9*sp - 9*pow(m_Lc,2) +
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))))/8.0;

                const double wm11z3 = (pow(m_Lb,-1)*pow(m_Lc,-1)*pow(-s0 + pow(m_Lb - m_Lc,2),4)*
                                      pow(8*pow(s0,3) - 8*pow(sp,3) +
                                      8*pow(sp,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      8*pow(s0,2)*(-3*sp +
                                      pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      sp*(-s0 + pow(m_Lb - m_Lc,2))*
                                      (3*s0 - 3*pow(m_Lb - m_Lc,2) +
                                      4*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      s0*(24*pow(sp,2) -
                                      16*sp*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      4*(-s0 + pow(m_Lb - m_Lc,2))*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      3*pow(-s0 + pow(m_Lb - m_Lc,2),2)) +
                                      6*(s0 - sp)*(2*s0 - 2*pow(m_Lb - m_Lc,2) +
                                      3*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))*
                                      pow(-s0 + pow(m_Lb - m_Lc,2),2)*
                                      pow(-s0 - 3*sp + 4*pow(m_Lb - m_Lc,2),-1),-1))/8.0;

                const double wm11z4 = (pow(m_Lb,-1)*pow(m_Lc,-1)*pow(s0 - sp,-2)*
                                      (30*m_Lc*pow(m_Lb,5) - 5*pow(m_Lb,6) - 5*pow(m_Lc,6) +
                                      64*sp*pow(s0,2) - 3*pow(s0,3) - 184*s0*pow(sp,2) +
                                      128*pow(sp,3) - 120*s0*sp*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      20*pow(s0,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      128*pow(sp,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      4*m_Lc*pow(m_Lb,3)*(65*s0 - 80*sp + 25*pow(m_Lc,2) -
                                      28*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      + pow(m_Lc,4)*(-65*s0 + 80*sp +
                                      28*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      + pow(m_Lb,4)*(-65*s0 + 80*sp - 75*pow(m_Lc,2) +
                                      28*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      + pow(m_Lc,2)*(-55*pow(s0,2) +
                                      80*s0*(3*sp + pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      8*sp*(25*sp + 17*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) +
                                      2*m_Lb*m_Lc*(15*pow(m_Lc,4) + 55*pow(s0,2) -
                                      80*s0*(3*sp + pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      8*sp*(25*sp + 17*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      2*pow(m_Lc,2)*(-65*s0 + 80*sp +
                                      28*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      )) + pow(m_Lb,2)*
                                      (-75*pow(m_Lc,4) - 55*pow(s0,2) +
                                      80*s0*(3*sp + pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      8*sp*(25*sp + 17*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      6*pow(m_Lc,2)*(-65*s0 + 80*sp +
                                      28*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)))))/32.0;

                const double wm11z5 = (pow(m_Lb,-1)*pow(m_Lc,-1)*pow(s0 - sp,-4)*
                                      (-((-s0 + sp)*(-112*s0*sp - 12*m_Lc*pow(m_Lb,3) +
                                      3*pow(m_Lb,4) + 6*(7*s0 - 8*sp)*pow(m_Lc,2) -
                                      12*m_Lb*m_Lc*(7*s0 - 8*sp + pow(m_Lc,2)) +
                                      6*pow(m_Lb,2)*(7*s0 - 8*sp + 3*pow(m_Lc,2)) +
                                      3*pow(m_Lc,4) + 35*pow(s0,2) + 80*pow(sp,2))*
                                      (-2*s0*pow(m_Lb - m_Lc,2) +
                                      2*sp*(s0 + pow(m_Lb - m_Lc,2)) - 2*pow(sp,2) -
                                      2*m_Lb*m_Lc*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lb,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lc,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) -
                                      2*(-s0 + pow(m_Lb - m_Lc,2))*
                                      (5*s0 - 8*sp + 3*pow(m_Lb - m_Lc,2))*pow(s0 - sp,2)*
                                      (-s0 + 3*pow(m_Lb - m_Lc,2) -
                                      2*(sp + pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) +
                                      2*(5*s0 - 6*sp + pow(m_Lb - m_Lc,2))*pow(s0 - sp,2)*
                                      pow(-s0 + pow(m_Lb - m_Lc,2),2) +
                                      2*pow(s0 - sp,2)*pow(-s0 + pow(m_Lb - m_Lc,2),3)))/32.0;

                const double wm12z1 = 0.0;

                const double wm12z2 = (pow(m_Lb,-2)*pow(m_Lc,-2)*pow(s0 - sp,-2)*
                                      pow(-2*s0*pow(m_Lb - m_Lc,2) +
                                      2*sp*(s0 + pow(m_Lb - m_Lc,2)) - 2*pow(sp,2) -
                                      2*m_Lb*m_Lc*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      s0*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lb,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lc,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5),2))/4.0;

                const double wm12z3 = (pow(m_Lb,-2)*pow(m_Lc,-2)*pow(s0 - sp,-2)*
                                      (-2*s0*pow(m_Lb - m_Lc,2) + 2*sp*(s0 + pow(m_Lb - m_Lc,2)) -
                                      2*pow(sp,2) - 2*m_Lb*m_Lc*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      s0*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lb,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lc,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))*
                                      (14*s0*sp + 12*m_Lc*pow(m_Lb,3) - 3*pow(m_Lb,4) -
                                      3*pow(m_Lc,4) - pow(s0,2) - 16*pow(sp,2) +
                                      6*s0*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      16*sp*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      4*m_Lb*m_Lc*(6*s0 - 9*sp + 3*pow(m_Lc,2) -
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      2*pow(m_Lc,2)*(-6*s0 + 9*sp +
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      2*pow(m_Lb,2)*(-6*s0 + 9*sp - 9*pow(m_Lc,2) +
                                      5*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))))/8.0;

                const double wm12z4 = (pow(m_Lb,-2)*pow(m_Lc,-2)*pow(s0 - sp,-2)*
                                      (-200*m_Lc*pow(m_Lb,7) + 25*pow(m_Lb,8) + 25*pow(m_Lc,8) -
                                      128*sp*pow(s0,3) + pow(s0,4) +
                                      1048*pow(s0,2)*pow(sp,2) - 2176*s0*pow(sp,3) +
                                      1280*pow(sp,4) + 440*sp*pow(s0,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      20*pow(s0,3)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) -
                                      1536*s0*pow(sp,2)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      1280*pow(sp,3)*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5) +
                                      pow(m_Lb,6)*(460*s0 - 560*sp + 700*pow(m_Lc,2) -
                                      164*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      - 4*pow(m_Lc,6)*(-115*s0 + 140*sp +
                                      41*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      + 8*m_Lc*pow(m_Lb,5)*
                                      (-345*s0 + 420*sp - 175*pow(m_Lc,2) +
                                      123*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5))\
                                      - 4*pow(m_Lc,2)*(-31*pow(s0,3) +
                                      32*pow(sp,2)*(23*sp +
                                      18*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + pow(s0,2)*
                                      (428*sp + 95*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      4*s0*sp*(277*sp +
                                      137*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) +
                                      2*pow(m_Lc,4)*(335*pow(s0,2) +
                                      4*sp*(275*sp + 151*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      2*s0*(680*sp + 179*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) +
                                      2*pow(m_Lb,4)*(875*pow(m_Lc,4) + 335*pow(s0,2) -
                                      30*pow(m_Lc,2)*(-115*s0 + 140*sp +
                                      41*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)) + 4*sp*(275*sp +
                                      151*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      2*s0*(680*sp + 179*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) -
                                      8*m_Lc*pow(m_Lb,3)*(175*pow(m_Lc,4) + 335*pow(s0,2) -
                                      10*pow(m_Lc,2)*(-115*s0 + 140*sp +
                                      41*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + 4*sp*(275*sp +
                                      51*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      2*s0*(680*sp + 179*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))) -
                                      4*pow(m_Lb,2)*(-175*pow(m_Lc,6) - 31*pow(s0,3) +
                                      32*pow(sp,2)*(23*sp +
                                      18*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + 15*pow(m_Lc,4)*
                                      (-115*s0 + 140*sp +
                                      41*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + pow(s0,2)*
                                      (428*sp + 95*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      4*s0*sp*(277*sp +
                                      137*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      3*pow(m_Lc,2)*(335*pow(s0,2) +
                                      4*sp*(275*sp +
                                      151*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      2*s0*(680*sp +
                                      179*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)))) +
                                      8*m_Lb*m_Lc*(-25*pow(m_Lc,6) - 31*pow(s0,3) +
                                      32*pow(sp,2)*(23*sp +
                                      18*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + 3*pow(m_Lc,4)*
                                      (-115*s0 + 140*sp +
                                      41*pow(-s0 + sp,0.5)*pow(sp - pow(m_Lb - m_Lc,2),0.5)
                                      ) + pow(s0,2)*
                                      (428*sp + 95*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) -
                                      4*s0*sp*(277*sp +
                                      137*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-335*pow(s0,2) -
                                      4*sp*(275*sp +
                                      151*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5)) +
                                      2*s0*(680*sp +
                                      179*pow(-s0 + sp,0.5)*
                                      pow(sp - pow(m_Lb - m_Lc,2),0.5))))))/64.0;

                const double wm12z5 = (pow(m_Lb,-2)*pow(m_Lc,-2)*pow(s0 - sp,-3)*
                                      (-93*sp*pow(s0,4) + pow(s0,5) + 908*pow(s0,3)*pow(sp,2) -
                                      2736*pow(s0,2)*pow(sp,3) + 3200*s0*pow(sp,4) -
                                      1280*pow(sp,5) + 340*sp*pow(s0,3)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      15*pow(s0,4)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      1616*pow(s0,2)*pow(sp,2)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      2560*s0*pow(sp,3)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      1280*pow(sp,4)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      pow(m_Lb,8)*(85*s0 - 85*sp -
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,8)*(85*s0 - 85*sp -
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      8*m_Lc*pow(m_Lb,7)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      4*pow(m_Lc,6)*(-365*s0*sp + 140*pow(s0,2) +
                                      225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) -
                                      4*pow(m_Lb,6)*(365*s0*sp - 140*pow(s0,2) -
                                      225*pow(sp,2) +
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      7*pow(m_Lc,2)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) + 4*pow(m_Lc,2)*
                                      (22*pow(s0,4) + 160*pow(sp,3)*
                                      (5*sp + 4*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - pow(s0,3)*(361*sp +
                                      70*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*s0*pow(sp,2)*
                                      (229*sp + 139*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + sp*pow(s0,2)*
                                      (1371*sp + 553*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) - 2*pow(m_Lc,4)*
                                      (-273*pow(s0,3) +
                                      7*pow(s0,2)*(237*sp +
                                      49*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*pow(sp,2)*
                                      (171*sp + 101*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 2*s0*sp*(1377*sp +
                                      559*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) + 8*pow(m_Lb,5)*
                                      (7*pow(m_Lc,3)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 3*m_Lc*(-365*s0*sp + 140*pow(s0,2) +
                                      225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)\
                                      + 83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) - 2*pow(m_Lb,4)*
                                      (-273*pow(s0,3) +
                                      35*pow(m_Lc,4)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 7*pow(s0,2)*
                                      (237*sp + 49*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*pow(sp,2)*
                                      (171*sp + 101*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 2*s0*sp*(1377*sp +
                                      559*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 30*pow(m_Lc,2)*
                                      (-365*s0*sp + 140*pow(s0,2) + 225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)\
                                      + 83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) + 8*m_Lc*pow(m_Lb,3)*
                                      (-273*pow(s0,3) +
                                      7*pow(m_Lc,4)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 7*pow(s0,2)*
                                      (237*sp + 49*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*pow(sp,2)*
                                      (171*sp + 101*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 2*s0*sp*(1377*sp +
                                      559*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 10*pow(m_Lc,2)*
                                      (-365*s0*sp + 140*pow(s0,2) + 225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)\
                                      + 83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) + 8*m_Lb*m_Lc*(-22*pow(s0,4) -
                                      160*pow(sp,3)*(5*sp +
                                      4*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + pow(m_Lc,6)*(-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + pow(s0,3)*(361*sp +
                                      70*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*s0*pow(sp,2)*
                                      (229*sp + 139*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - sp*pow(s0,2)*
                                      (1371*sp + 553*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 3*pow(m_Lc,4)*
                                      (-365*s0*sp + 140*pow(s0,2) + 225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)\
                                      + 83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + pow(m_Lc,2)*(-273*pow(s0,3) +
                                      7*pow(s0,2)*(237*sp +
                                      49*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)) +
                                      8*pow(sp,2)*(171*sp +
                                      101*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)) -
                                      2*s0*sp*(1377*sp +
                                      559*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)))) -
                                      4*pow(m_Lb,2)*(-22*pow(s0,4) -
                                      160*pow(sp,3)*(5*sp +
                                      4*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 7*pow(m_Lc,6)*
                                      (-85*s0 + 85*sp +
                                      11*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + pow(s0,3)*(361*sp +
                                      70*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*s0*pow(sp,2)*
                                      (229*sp + 139*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - sp*pow(s0,2)*
                                      (1371*sp + 553*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 15*pow(m_Lc,4)*
                                      (-365*s0*sp + 140*pow(s0,2) + 225*pow(sp,2) -
                                      72*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)\
                                      + 83*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 3*pow(m_Lc,2)*
                                      (-273*pow(s0,3) +
                                      7*pow(s0,2)*(237*sp +
                                      49*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)) +
                                      8*pow(sp,2)*(171*sp +
                                      101*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)) -
                                      2*s0*sp*(1377*sp +
                                      559*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5))))))/32.0;

                const double wm13z1 = 0.0;

                const double wm13z2 = 0.0;

                const double wm13z3 = (pow(m_Lb,-3)*pow(m_Lc,-3)*pow(s0 - sp,-3)*
                                      pow(2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)),3))/8.0;

                const double wm13z4 = (3*pow(m_Lb,-3)*pow(m_Lc,-3)*pow(s0 - sp,-3)*
                                      (14*s0*sp + 12*m_Lc*pow(m_Lb,3) - 3*pow(m_Lb,4) -
                                      3*pow(m_Lc,4) - pow(s0,2) - 16*pow(sp,2) +
                                      6*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      16*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      4*m_Lb*m_Lc*(6*s0 - 9*sp + 3*pow(m_Lc,2) -
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      2*pow(m_Lc,2)*(-6*s0 + 9*sp +
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      2*pow(m_Lb,2)*(-6*s0 + 9*sp - 9*pow(m_Lc,2) +
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)))*
                                      pow(2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)),2))/32.0;

                const double wm13z5 = (-3*pow(m_Lb,-3)*pow(m_Lc,-3)*pow(s0 - sp,-3)*
                                      (2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)))*
                                      (136*m_Lc*pow(m_Lb,7) - 17*pow(m_Lb,8) - 17*pow(m_Lc,8) +
                                      96*sp*pow(s0,3) - pow(s0,4) - 752*pow(s0,2)*pow(sp,2) +
                                      1536*s0*pow(sp,3) - 896*pow(sp,4) -
                                      320*sp*pow(s0,2)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      16*pow(s0,3)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      1088*s0*pow(sp,2)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      896*pow(sp,3)*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      8*m_Lc*pow(m_Lb,5)*(237*s0 - 288*sp + 119*pow(m_Lc,2) -
                                      84*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      4*pow(m_Lc,6)*(-79*s0 + 96*sp +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      4*pow(m_Lb,6)*(-79*s0 + 96*sp - 119*pow(m_Lc,2) +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      4*pow(m_Lc,2)*(-23*pow(s0,3) +
                                      4*pow(s0,2)*(76*sp +
                                      17*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 16*pow(sp,2)*
                                      (32*sp + 25*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*s0*sp*(97*sp +
                                      48*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) - 8*m_Lc*pow(m_Lb,3)*
                                      (-119*pow(m_Lc,4) - 235*pow(s0,2) +
                                      10*pow(m_Lc,2)*(-79*s0 + 96*sp +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*s0*(118*sp +
                                      31*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*sp*(95*sp +
                                      52*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) + 2*pow(m_Lb,4)*
                                      (-595*pow(m_Lc,4) - 235*pow(s0,2) +
                                      30*pow(m_Lc,2)*(-79*s0 + 96*sp +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*s0*(118*sp +
                                      31*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*sp*(95*sp +
                                      52*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) - 2*pow(m_Lc,4)*
                                      (235*pow(s0,2) - 8*s0*
                                      (118*sp + 31*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 8*sp*(95*sp +
                                      52*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      ) - 8*m_Lb*m_Lc*(-17*pow(m_Lc,6) - 23*pow(s0,3) +
                                      4*pow(s0,2)*(76*sp +
                                      17*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 16*pow(sp,2)*
                                      (32*sp + 25*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 3*pow(m_Lc,4)*
                                      (-79*s0 + 96*sp +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*s0*sp*(97*sp +
                                      48*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + pow(m_Lc,2)*(-235*pow(s0,2) +
                                      8*s0*(118*sp +
                                      31*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)) -
                                      8*sp*(95*sp +
                                      52*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5)))) +
                                      4*pow(m_Lb,2)*(-119*pow(m_Lc,6) - 23*pow(s0,3) +
                                      4*pow(s0,2)*(76*sp +
                                      17*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 16*pow(sp,2)*
                                      (32*sp + 25*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      + 15*pow(m_Lc,4)*
                                      (-79*s0 + 96*sp +
                                      28*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 8*s0*sp*(97*sp +
                                      48*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5))
                                      - 3*pow(m_Lc,2)*
                                      (235*pow(s0,2) -
                                      8*s0*(118*sp +
                                      31*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      8*sp*(95*sp +
                                      52*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),
                                      0.5))))))/128.0;

                const double wm14z1 = 0.0;

                const double wm14z2 = 0.0;

                const double wm14z3 = 0.0;

                const double wm14z4 = (pow(m_Lb,-4)*pow(m_Lc,-4)*pow(s0 - sp,-4)*
                                      pow(2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)),4))/16.0;

                const double wm14z5 = (pow(m_Lb,-4)*pow(m_Lc,-4)*pow(s0 - sp,-4)*
                                      (14*s0*sp + 12*m_Lc*pow(m_Lb,3) - 3*pow(m_Lb,4) -
                                      3*pow(m_Lc,4) - pow(s0,2) - 16*pow(sp,2) +
                                      6*s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      16*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) +
                                      4*m_Lb*m_Lc*(6*s0 - 9*sp + 3*pow(m_Lc,2) -
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      2*pow(m_Lc,2)*(-6*s0 + 9*sp +
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      2*pow(m_Lb,2)*(-6*s0 + 9*sp - 9*pow(m_Lc,2) +
                                      5*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)))*
                                      pow(2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)),3))/16.0;

                const double wm15z1 = 0.0;

                const double wm15z2 = 0.0;

                const double wm15z3 = 0.0;

                const double wm15z4 = 0.0;

                const double wm15z5 = (pow(m_Lb,-5)*pow(m_Lc,-5)*pow(s0 - sp,-5)*
                                      pow(2*s0*sp - 2*pow(sp,2) +
                                      s0*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*sp*pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5) -
                                      2*m_Lb*m_Lc*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lb,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)) +
                                      pow(m_Lc,2)*(-2*s0 + 2*sp +
                                      pow(-s0 + sp,0.5)*
                                      pow(2*m_Lb*m_Lc + sp - pow(m_Lb,2) - pow(m_Lc,2),0.5)),5))/32.0;

                const double wm11 =  wm11z1 * z + wm11z2 * z2 + wm11z3 * z3 + wm11z4 * z4 + wm11z5 * z5;

                const double wm12 =  wm12z1 * z + wm12z2 * z2 + wm12z3 * z3 + wm12z4 * z4 + wm12z5 * z5;

                const double wm13 =  wm13z1 * z + wm13z2 * z2 + wm13z3 * z3 + wm13z4 * z4 + wm13z5 * z5;

                const double wm14 =  wm14z1 * z + wm14z2 * z2 + wm14z3 * z3 + wm14z4 * z4 + wm14z5 * z5;

                const double wm15 =  wm15z1 * z + wm15z2 * z2 + wm15z3 * z3 + wm15z4 * z4 + wm15z5 * z5;

                return 1.0
                    + _zetapone             * wm11
                    + _zetappone    / 2.0   * wm12
                    + _zetapppone   / 6.0   * wm13
                    + _zetappppone  / 24.0  * wm14
                    + _zetapppppone / 120.0 * wm15;
            }

            Diagnostics diagnostics() const
            {
                Diagnostics results;

                // Inputs
                {
                    results.add(Diagnostics::Entry{ _q2(1.0), "q2(w = 1.0)" });
                    results.add(Diagnostics::Entry{ _zeta_power_series(1.0), "q2(w = 1.0)" });
                }

                return results;
            }
    };
}


#endif