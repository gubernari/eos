/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2020 Danny van Dyk
 * Copyright (c) 2020 Nico Gubernari
 * Copyright (c) 2020 Christoph Bobeth
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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_IMPL_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_IMPL_HH 1

#include <eos/form-factors/parametric-bgl1997.hh>
#include <eos/utils/kinematic.hh>
#include <eos/utils/model.hh>
#include <eos/utils/power_of.hh>

#include <gsl/gsl_sf_dilog.h>

namespace eos
{
    double
    BGL1997FormFactors<BToDstar>::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

    double
    BGL1997FormFactors<BToDstar>::_chi_T(const double & u) const
    {
        // [BGL1997] eq.(4.1) + (4.2) + (4.5, 4.8)
        // limit for u -> 1 in eq.(4.11)
        // limit for u -> 0 in eq.(4.10)
        const double aS = _model->alpha_s(_mu) / M_PI;
        const double u2 = power_of<2>(u);
        const double u3 = power_of<3>(u);
        const double u4 = power_of<4>(u);
        const double u5 = power_of<5>(u);
        const double u6 = power_of<6>(u);
        const double lnu = std::log(u2);
       
        const double chi_pert =
            1.0 / 32.0 * ( 
               (1.0 - u2) * (3. + 4. * u - 21. * u2 + 40. * u3 - 21. * u4 + 4. * u5 + 3. * u6)
               + 12. * u3 * (2. - 3. * u + 2. * u2) * lnu 
            )
            + aS / (576.0 * (1.0 - u2)) * (
               power_of<2>(1.0 - u2) * (75. + 360. * u - 1031. * u2 + 1776. * u3 - 1031. * u4 + 360. * u5 + 75. * u6)
               + 4.0 * u * (1.0 - u2) * (18. - 99. * u + 732. * u2 - 1010. * u3 + 732. * u4 - 99. * u5 + 18. * u6) * lnu
               + 4.0 * u3 * (108. - 324. * u + 648. * u2 - 456. * u3 + 132. * u4 + 59. * u5 - 12. * u6 - 9. *u6 * u) * lnu * lnu
               + 8.0 * power_of<3>(1.0 - u2) * (9. + 12. * u - 32. * u2 + 12. * u3 + 9. * u4) * gsl_sf_dilog(1.0 - u2)
            );
       
        const double chi_cond =
            -1.0 * _cond_qq / 2.0 * (2. - 3. * u + 2. * u2) 
            -1.0 * _cond_G2 / (24.0 * _m_b_pole() * power_of<2>(1.0 - u2) ) * (
                (1.0 - u2) * (2. - 104. * u + 148. * u2 - 270. * u3 + 145. * u4 - 104. * u5 + 5. * u6 - 2. * u6 * u)
               - 12.0 * u * (3. - 5. * u + 17. * u2 - 15. * u3 + 17. * u4 - 5. * u5 + 3. * u6) * lnu
            );
        
        //  TODO catch u->1 and u->0, but should not be needed for b->c
       
        return chi_pert / (power_of<2>(_m_b_pole() * M_PI) * power_of<5>(1.0 - u2))
             + chi_cond / (power_of<5>(_m_b_pole() * (1.0 - u2)) );
    }

    double
    BGL1997FormFactors<BToDstar>::_chi_L(const double & u) const
    {
        // [BGL1997] eq.(4.1) + (4.3) + (4.6, 4.9)
        // limit for u -> 1 in eq.(4.11)
        // limit for u -> 0 in eq.(4.10)
        const double aS = _model->alpha_s(_mu) / M_PI;
        const double u2 = power_of<2>(u);
        const double u3 = power_of<3>(u);
        const double u4 = power_of<4>(u);
        const double u5 = power_of<5>(u);
        const double lnu = std::log(u2);
       
        const double chi_pert =
            1.0 / 8.0 * ( 
               (1.0 - u2) * (1. + u + u2) * (1 - 4. * u + u2) - 6. * u3 * lnu 
            )
            + aS / (48.0 * (1.0 - u2)) * (
               power_of<2>(1.0 - u2) * (1. - 36. * u - 22. * u2 - 36. * u3 + u4)
               - 2.0 * u * (1.0 - u2) * (9. + 4. * u + 66. * u2 + 4. * u3 + 9. * u4) * lnu
               - 4.0 * u3 * (9. + 18. * u2 - 2. * u3 - 3. * u4 + u5) * lnu * lnu
               + 8.0 * power_of<3>(1.0 - u2) * (1. - 3. * u + u2) * gsl_sf_dilog(1.0 - u2)
            );
       
        const double chi_cond =
            _cond_qq + _cond_G2 / (12.0 * _m_b_pole() * power_of<2>(1.0 - u2) ) * (
                (1.0 - u2) * (1. - 21. * u + 10. * u2 - 20. * u3 + u4 - u5)
               - 3.0 * u * (3. - 2. * u + 8. * u2 - 2. * u3 + 3. * u4) * lnu
            );

        //  TODO catch u->1 and u->0, but should not be needed for b->c
       
        return chi_pert/(M_PI * M_PI * power_of<3>(1.0 - u2))
             + chi_cond / (power_of<3>(_m_b_pole() * (1.0 - u2)) );    
    }

    double
    BGL1997FormFactors<BToDstar>::_phi(const double & s, const unsigned & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const
    { 
        const double sq_tp    = std::sqrt(_t_p);
        const double sq_tp_t  = std::sqrt(_t_p - s);
        const double sq_tp_t0 = std::sqrt(_t_p - _t_0);
        const double sq_tp_tm = std::sqrt(_t_p - _t_m);

        // [BGL1997] eq. (4.14) for OPE at Q^2 = -q^2 = 0 
        // => generalization for q^2 != 0 possible, see eq.(4.15)      
        return std::sqrt(nI / (K * M_PI * chi)) * (sq_tp_t + sq_tp_t0)
               * std::pow(sq_tp_t / sq_tp_t0, 1.0 / 4.0)
               * std::pow(_t_p - s, a / 4.0)
               * std::pow(sq_tp_t + sq_tp_tm, b / 2.0)
               * 1.0 / std::pow(sq_tp_t + sq_tp, c + 3.0);
    }

    std::string
    BGL1997FormFactors<BToDstar>::_par_name(const std::string & ff_name)
    {
        return std::string("B->D^*") + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    BGL1997FormFactors<BToDstar>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        _model(Model::make("SM", p, o)), 
        _a_g{{   UsedParameter(p[_par_name("g_0")],  *this),
                 UsedParameter(p[_par_name("g_1")],  *this),
                 UsedParameter(p[_par_name("g_2")],  *this),
                 UsedParameter(p[_par_name("g_3")],  *this) }},
        _a_f{{   UsedParameter(p[_par_name("f_0")],  *this),
                 UsedParameter(p[_par_name("f_1")],  *this),
                 UsedParameter(p[_par_name("f_2")],  *this),
                 UsedParameter(p[_par_name("f_3")],  *this) }},
        _a_F1{{  UsedParameter(p[_par_name("F1_0")], *this),
                 UsedParameter(p[_par_name("F1_1")], *this),
                 UsedParameter(p[_par_name("F1_2")], *this),
                 UsedParameter(p[_par_name("F1_3")], *this) }},
        _a_F2{{  UsedParameter(p[_par_name("F2_0")], *this),
                 UsedParameter(p[_par_name("F2_1")], *this),
                 UsedParameter(p[_par_name("F2_2")], *this),
                 UsedParameter(p[_par_name("F2_3")], *this) }},
        _t_0(p["B->D^*::t_0@BGL1997"], *this),
        _m_c_pole(p["B->D^*::m_c^pole@BGL1997"], *this),
        _m_b_pole(p["B->D^*::m_b^pole@BGL1997"], *this),
        _cond_qq(p["B->D^*::<qq>@BGL1997"], *this),         // [BGL1997] quark condensate
        _cond_G2(p["B->D^*::<alS/pi G^2>@BGL1997"], *this), // [BGL1997] gluon condensate
        _mu(4.2),                                           // TODO remove hard-coded numerical value
        _mB(BToDstar::mB),
        _mB2(power_of<2>(_mB)),
        _mV(BToDstar::mV),
        _mV2(power_of<2>(_mV)),
        _t_p(power_of<2>(_mB + _mV)),
        _t_m(power_of<2>(_mB - _mV)),
        _chi_T_pl(_chi_T(_m_c_hat())),  // TODO probably wrong, because depend on UsedParameter's and need to be updated during analysis
        _chi_L_pl(_chi_L(_m_c_hat())),
        _chi_T_mi(_chi_T(-1.0 * _m_c_hat())),
        _chi_L_mi(_chi_L(-1.0 * _m_c_hat()))
    {
    }

    BGL1997FormFactors<BToDstar>::~BGL1997FormFactors() = default;

    FormFactors<PToV> *
    BGL1997FormFactors<BToDstar>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    double
    BGL1997FormFactors<BToDstar>::g(const double & s) const
    {
        // resonances for 1^- at 6.329, 6.910, 7.020 from [BGS2017] table III with change 6.920->6.910
        const double blaschke = _z(s, 6.239) * _z(s, 6.910) * _z(s, 7.020);
        const double phi      = _phi(s, 96, 3, 3, 1, _chi_T_pl);
        const double z        = _z(s, _t_0);
        const double series   = _a_g[0] + _a_g[1] * z + _a_g[2] * z * z + _a_g[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::f(const double & s) const
    {
        // resonances for 1^+ at 6.739, 6.750, 7.145, 7.150 
        const double blaschke = _z(s, 6.739) * _z(s, 6.750) * _z(s, 7.145) * _z(s, 7.150);
        const double phi      = _phi(s, 24, 1, 1, 1, _chi_T_mi);
        const double z        = _z(s, _t_0);
        const double series   = _a_f[0] + _a_f[1] * z + _a_f[2] * z * z + _a_f[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::F1(const double & s) const
    {
        // resonances for 1^+ at 6.739, 6.750, 7.145, 7.150 
        const double blaschke = _z(s, 6.739) * _z(s, 6.750) * _z(s, 7.145) * _z(s, 7.150);  //  TODO
        const double phi      = _phi(s, 48, 1, 1, 2, _chi_T_mi);
        const double z        = _z(s, _t_0);
        const double series   = _a_F1[0] + _a_F1[1] * z + _a_F1[2] * z * z + _a_F1[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::F2(const double & s) const
    {
        // resonances for 0^- at 6.275, 6.842, 7.250
        const double blaschke = _z(s, 6.275) * _z(s, 6.842) * _z(s, 7.250); //  TODO
        const double phi      = _phi(s, 64, 3, 3, 1, _chi_L_mi);
        const double z        = _z(s, _t_0);
        const double series   = _a_F2[0] + _a_F2[1] * z + _a_F2[2] * z * z + _a_F2[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::v(const double & s) const
    {
        return (_mB + _mV) / 2.0 * g(s);
    }

    double
    BGL1997FormFactors<BToDstar>::a_0(const double & s) const
    {
        return F2(s) / 2.0;
    }

    double
    BGL1997FormFactors<BToDstar>::a_1(const double & s) const
    {
        return -1.0/(_mB + _mV) * f(s);
    }

    double
    BGL1997FormFactors<BToDstar>::a_2(const double & s) const
    {
        return (_mB + _mV) / eos::lambda(_mB2, _mV2, s) * ((_mB2 - _mV2 - s) * f(s) - 2.0 * _mV * F1(s));
    }

    double
    BGL1997FormFactors<BToDstar>::a_12(const double & s) const
    {
        return F1(s) / (8.0 * _mB * _mV);
    }

    double
    BGL1997FormFactors<BToDstar>::t_1(const double & s) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_2(const double & s) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_3(const double & s) const
    {
        return 0.0;  //  TODO
    }

    double
    BGL1997FormFactors<BToDstar>::t_23(const double & s) const
    {
        return 0.0;  //  TODO
    }


   
   



// TODO same function used in <BToDstar>
    double
    BGL1997FormFactors<BToD>::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

// TODO same function used in <BToDstar>
    double
    BGL1997FormFactors<BToD>::_chi_T(const double & u) const
    {
        // [BGL1997] eq.(4.1) + (4.2) + (4.5, 4.8)
        // limit for u -> 1 in eq.(4.11)
        // limit for u -> 0 in eq.(4.10)
        const double aS = _model->alpha_s(_mu) / M_PI;
        const double u2 = power_of<2>(u);
        const double u3 = power_of<3>(u);
        const double u4 = power_of<4>(u);
        const double u5 = power_of<5>(u);
        const double u6 = power_of<6>(u);
        const double lnu = std::log(u2);
       
        const double chi_pert =
            1.0 / 32.0 * ( 
               (1.0 - u2) * (3. + 4. * u - 21. * u2 + 40. * u3 - 21. * u4 + 4. * u5 + 3. * u6)
               + 12. * u3 * (2. - 3. * u + 2. * u2) * lnu 
            )
            + aS / (576.0 * (1.0 - u2)) * (
               power_of<2>(1.0 - u2) * (75. + 360. * u - 1031. * u2 + 1776. * u3 - 1031. * u4 + 360. * u5 + 75. * u6)
               + 4.0 * u * (1.0 - u2) * (18. - 99. * u + 732. * u2 - 1010. * u3 + 732. * u4 - 99. * u5 + 18. * u6) * lnu
               + 4.0 * u3 * (108. - 324. * u + 648. * u2 - 456. * u3 + 132. * u4 + 59. * u5 - 12. * u6 - 9. *u6 * u) * lnu * lnu
               + 8.0 * power_of<3>(1.0 - u2) * (9. + 12. * u - 32. * u2 + 12. * u3 + 9. * u4) * gsl_sf_dilog(1.0 - u2)
            );
       
        const double chi_cond =
            -1.0 * _cond_qq / 2.0 * (2. - 3. * u + 2. * u2) 
            -1.0 * _cond_G2 / (24.0 * _m_b_pole() * power_of<2>(1.0 - u2) ) * (
                (1.0 - u2) * (2. - 104. * u + 148. * u2 - 270. * u3 + 145. * u4 - 104. * u5 + 5. * u6 - 2. * u6 * u)
               - 12.0 * u * (3. - 5. * u + 17. * u2 - 15. * u3 + 17. * u4 - 5. * u5 + 3. * u6) * lnu
            );
        
        //  TODO catch u->1 and u->0, but should not be needed for b->c
       
        return chi_pert / (power_of<2>(_m_b_pole() * M_PI) * power_of<5>(1.0 - u2))
             + chi_cond / (power_of<5>(_m_b_pole() * (1.0 - u2)) );
    }

// TODO same function used in <BToDstar>
    double
    BGL1997FormFactors<BToD>::_chi_L(const double & u) const
    {
        // [BGL1997] eq.(4.1) + (4.3) + (4.6, 4.9)
        // limit for u -> 1 in eq.(4.11)
        // limit for u -> 0 in eq.(4.10)
        const double aS = _model->alpha_s(_mu) / M_PI;
        const double u2 = power_of<2>(u);
        const double u3 = power_of<3>(u);
        const double u4 = power_of<4>(u);
        const double u5 = power_of<5>(u);
        const double lnu = std::log(u2);
       
        const double chi_pert =
            1.0 / 8.0 * ( 
               (1.0 - u2) * (1. + u + u2) * (1 - 4. * u + u2) - 6. * u3 * lnu 
            )
            + aS / (48.0 * (1.0 - u2)) * (
               power_of<2>(1.0 - u2) * (1. - 36. * u - 22. * u2 - 36. * u3 + u4)
               - 2.0 * u * (1.0 - u2) * (9. + 4. * u + 66. * u2 + 4. * u3 + 9. * u4) * lnu
               - 4.0 * u3 * (9. + 18. * u2 - 2. * u3 - 3. * u4 + u5) * lnu * lnu
               + 8.0 * power_of<3>(1.0 - u2) * (1. - 3. * u + u2) * gsl_sf_dilog(1.0 - u2)
            );
       
        const double chi_cond =
            _cond_qq + _cond_G2 / (12.0 * _m_b_pole() * power_of<2>(1.0 - u2) ) * (
                (1.0 - u2) * (1. - 21. * u + 10. * u2 - 20. * u3 + u4 - u5)
               - 3.0 * u * (3. - 2. * u + 8. * u2 - 2. * u3 + 3. * u4) * lnu
            );

        //  TODO catch u->1 and u->0, but should not be needed for b->c
       
        return chi_pert/(M_PI * M_PI * power_of<3>(1.0 - u2))
             + chi_cond / (power_of<3>(_m_b_pole() * (1.0 - u2)) );    
    }

// TODO same function used in <BToDstar>
    double
    BGL1997FormFactors<BToD>::_phi(const double & s, const unsigned & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const
    { 
        const double sq_tp    = std::sqrt(_t_p);
        const double sq_tp_t  = std::sqrt(_t_p - s);
        const double sq_tp_t0 = std::sqrt(_t_p - _t_0);
        const double sq_tp_tm = std::sqrt(_t_p - _t_m);

        // [BGL1997] eq. (4.14) for OPE at Q^2 = -q^2 = 0 
        // => generalization for q^2 != 0 possible, see eq.(4.15)      
        return std::sqrt(nI / (K * M_PI * chi)) * (sq_tp_t + sq_tp_t0)
               * std::pow(sq_tp_t / sq_tp_t0, 1.0 / 4.0)
               * std::pow(_t_p - s, a / 4.0)
               * std::pow(sq_tp_t + sq_tp_tm, b / 2.0)
               * 1.0 / std::pow(sq_tp_t + sq_tp, c + 3.0);
    }

    std::string
    BGL1997FormFactors<BToD>::_par_name(const std::string & ff_name)
    {
        return std::string("B->D") + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    BGL1997FormFactors<BToD>::BGL1997FormFactors(const Parameters & p, const Options & o) :
        _model(Model::make("SM", p, o)), 
        _a_f_p{{ UsedParameter(p[_par_name("f_p_0")], *this),
                 UsedParameter(p[_par_name("f_p_1")], *this),
                 UsedParameter(p[_par_name("f_p_2")], *this),
                 UsedParameter(p[_par_name("f_p_3")], *this) }},
        _a_f_0{{ UsedParameter(p[_par_name("f_0_0")], *this),
                 UsedParameter(p[_par_name("f_0_1")], *this),
                 UsedParameter(p[_par_name("f_0_2")], *this),
                 UsedParameter(p[_par_name("f_0_3")], *this) }},
        _a_f_t{{ UsedParameter(p[_par_name("f_t_0")], *this),
                 UsedParameter(p[_par_name("f_t_1")], *this),
                 UsedParameter(p[_par_name("f_t_2")], *this),
                 UsedParameter(p[_par_name("f_t_3")], *this) }},
        _t_0(p["B->D::t_0@BGL1997"], *this),
        _m_c_pole(p["B->D::m_c^pole@BGL1997"], *this),    // TODO <- same as in <BToDstar>
        _m_b_pole(p["B->D::m_b^pole@BGL1997"], *this),    // TODO <- same as in <BToDstar>
        _cond_qq(p["B->D::<qq>@BGL1997"], *this),         // [BGL1997] quark condensate  TODO <- same as in <BToDstar>
        _cond_G2(p["B->D::<alS/pi G^2>@BGL1997"], *this), // [BGL1997] gluon condensate  TODO <- same as in <BToDstar>
        _mu(4.2),                                         // TODO remove hard-coded numerical value  TODO <- same as in <BToDstar>
        _mB(BToD::m_B),
        _mB2(power_of<2>(_mB)),
        _mP(BToD::m_P),
        _mP2(power_of<2>(_mP)),
        _t_p(power_of<2>(_mB + _mP)),
        _t_m(power_of<2>(_mB - _mP)),
        _chi_T_pl(_chi_T(_m_c_hat())),  // TODO probably wrong, because depend on UsedParameter's and need to be updated during analysis
        _chi_L_pl(_chi_L(_m_c_hat())),
        _chi_T_mi(_chi_T(-1.0 * _m_c_hat())),
        _chi_L_mi(_chi_L(-1.0 * _m_c_hat()))
    {
    }

    BGL1997FormFactors<BToD>::~BGL1997FormFactors() = default;

    FormFactors<PToP> *
    BGL1997FormFactors<BToD>::make(const Parameters & parameters, const Options & options)
    {
        return new BGL1997FormFactors(parameters, options);
    }

    double
    BGL1997FormFactors<BToD>::f_p(const double & s) const
    {
        // resonances for 1^- at 6.329, 6.910, 7.020 from [BGS2017] table III with change 6.920->6.910
        const double blaschke = _z(s, 6.239) * _z(s, 6.910) * _z(s, 7.020);
        const double phi      = _phi(s, 48, 3, 3, 2, _chi_T_pl);
        const double z        = _z(s, _t_0);
        const double series   = _a_f_p[0] + _a_f_p[1] * z + _a_f_p[2] * z * z + _a_f_p[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToD>::f_0(const double & s) const
    {
        // resonances for 0^+ at 6.704, 7.122 
        const double blaschke = _z(s, 6.704) * _z(s, 7.122) * _z(s, 7.145) * _z(s, 7.150);
        const double phi      = _phi(s, 16, 1, 1, 1, _chi_L_pl);
        const double z        = _z(s, _t_0);
        const double series   = _a_f_0[0] + _a_f_0[1] * z + _a_f_0[2] * z * z + _a_f_0[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToD>::f_t(const double & s) const
    {
        return 0.0; //  TODO
    }

}

#endif
