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
#include <eos/utils/power_of.hh>

namespace eos
{
    double
    BGL1997FormFactors<BToDstar>::_z(const double & t, const double & t_0) const
    {
        return (std::sqrt(_t_p - t) - std::sqrt(_t_p - t_0)) / (std::sqrt(_t_p - t) + std::sqrt(_t_p - t_0));
    }

    double
    BGL1997FormFactors<BToDstar>::_phi(const double & s, const double & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const
    {
        // eq. (4.14)
        return 0.0;
    }

    std::string
    BGL1997FormFactors<BToDstar>::_par_name(const std::string & ff_name)
    {
        return std::string("B->D^*") + std::string("::a^") + ff_name + std::string("@BGL1997");
    }

    BGL1997FormFactors<BToDstar>::BGL1997FormFactors(const Parameters & p, const Options & /*o*/) :
        _a_g{{   UsedParameter(p[_par_name("g_0")],  *this),
                 UsedParameter(p[_par_name("g_1")],  *this),
                 UsedParameter(p[_par_name("g_2")],  *this),
                 UsedParameter(p[_par_name("g_3")],  *this)}},
        _a_f{{   UsedParameter(p[_par_name("f_0")],  *this),
                 UsedParameter(p[_par_name("f_1")],  *this),
                 UsedParameter(p[_par_name("f_2")],  *this),
                 UsedParameter(p[_par_name("f_3")],  *this) }},
        _a_F1{{  UsedParameter(p[_par_name("F1_0")],   *this),
                 UsedParameter(p[_par_name("F1_1")],   *this),
                 UsedParameter(p[_par_name("F1_2")],   *this),
                 UsedParameter(p[_par_name("F1_3")],   *this) }},
        _a_F2{{  UsedParameter(p[_par_name("F2_0")],  *this),
                 UsedParameter(p[_par_name("F2_1")],  *this),
                 UsedParameter(p[_par_name("F2_2")],  *this),
                 UsedParameter(p[_par_name("F2_3")],  *this) }},
        _t_0(p["B->D^*::t_0@BGL1997"], *this),
        _mB(BToDstar::mB),
        _mB2(power_of<2>(_mB)),
        _mV(BToDstar::mV),
        _mV2(power_of<2>(_mV)),
        _t_p(power_of<2>(_mB + _mV)),
        _t_m(power_of<2>(_mB - _mV))
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
        // resonances at 6.329, 6.910, 7.020 
        const double blaschke = _z(s, 6.239) * _z(s, 6.910) * _z(s, 7.020);
        const double phi      = _phi(s, 96, 3, 3, 1, chi_1m);
        const double z        = _z(s, _t_0);
        const double series   = _a_g[0] + _a_g[1] * z + _a_g[2] * z * z + _a_g[3] * z * z * z;

        return series / phi / blaschke;
    }

    double
    BGL1997FormFactors<BToDstar>::f(const double & s) const
    {
        // resonances at 6.739, 6.750, 7.145, 7.150 
        const double blaschke = _z(s, 6.739) * _z(s, 6.750) * _z(s, 7.145)  * _z(s, 7.150);
        const double phi      = _phi(s, 96, 3, 3, 1, chi_1m);
        const double z        = _z(s, _t_0);
        const double series   = _a_f[0] + _a_f[1] * z + _a_f[2] * z * z + _a_f[3] * z * z * z;

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
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::a_1(const double & s) const
    {
        return -1.0/(_mB + _mV) * f(s);
    }

    double
    BGL1997FormFactors<BToDstar>::a_2(const double & s) const
    {
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::a_12(const double & s) const
    {
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::t_1(const double & s) const
    {
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::t_2(const double & s) const
    {
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::t_3(const double & s) const
    {
        return 0.0;
    }

    double
    BGL1997FormFactors<BToDstar>::t_23(const double & s) const
    {
        return 0.0;
    }
}

#endif
