/* vim: set sw=4 sts=4 et tw=120 foldmethod=syntax : */

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

#ifndef EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH
#define EOS_GUARD_EOS_FORM_FACTORS_PARAMETRIC_BGL1997_HH 1

#include <eos/form-factors/form-factors-fwd.hh>
#include <eos/form-factors/mesonic.hh>
#include <eos/form-factors/mesonic-processes.hh>
#include <eos/utils/model.hh>

// TODO move references
// [BGL1997] C.G. Boyd, B. Grinstein and R.F. Lebed, hep-ph/9705252
// [BGS2017] D. Bigi, P. Gambino and S. Schacht, arXiv:1707.09509

namespace eos
{
    template <typename Process_> class BGL1997FormFactors;

    template <>
    class BGL1997FormFactors<BToDstar> :
        public FormFactors<PToV>
    {
        private:
            
            std::shared_ptr<Model> _model;
            
            std::array<UsedParameter, 4> _a_g, _a_f, _a_F1, _a_F2;
            UsedParameter _t_0;
            UsedParameter _m_c_pole, _m_b_pole;
            UsedParameter _cond_qq, _cond_G2;

            const double _mu;   // TODO maybe declare the scale of alpha_s(mu) as a UsedParameter
            const double _mB, _mB2, _mV, _mV2;
            const double _t_p, _t_m;
            const double _chi_T_pl, _chi_L_pl;
            const double _chi_T_mi, _chi_L_mi;

            static constexpr double chi_1m = 1.0;
            static constexpr double chi_1p = 1.0;
            static constexpr double nI = 2.0;     // [BGL1997] isospin Clebsch-Gordan factor for B->D* and B->D

            double inline _m_c_hat() const { return _m_c_pole() / _m_b_pole(); }

            double _z(const double & t, const double & t_0) const;
            
            double _chi_T(const double & u) const;

            double _chi_L(const double & u) const;

            double _phi(const double & s, const unsigned & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const;

            static std::string _par_name(const std::string & ff_name);

        public:

            BGL1997FormFactors(const Parameters &, const Options &);

            ~BGL1997FormFactors();

            static FormFactors<PToV> * make(const Parameters & parameters, const Options & options);

            double g(const double & s) const;

            double f(const double & s) const;

            double F1(const double & s) const;

            double F2(const double & s) const;

            virtual double v(const double & s) const;

            virtual double a_0(const double & s) const;

            virtual double a_1(const double & s) const;

            virtual double a_2(const double & s) const;

            virtual double a_12(const double & s) const;

            virtual double t_1(const double & s) const;

            virtual double t_2(const double & s) const;

            virtual double t_3(const double & s) const;

            virtual double t_23(const double & s) const;
    };



    template <>
    class BGL1997FormFactors<BToD> :
        public FormFactors<PToP>
    {
        private:
            
            std::shared_ptr<Model> _model;
            
            std::array<UsedParameter, 4> _a_f_p, _a_f_0, _a_f_t;
            UsedParameter _t_0;
            UsedParameter _m_c_pole, _m_b_pole;
            UsedParameter _cond_qq, _cond_G2;

            const double _mu;   // TODO maybe declare the scale of alpha_s(mu) as a UsedParameter
            const double _mB, _mB2, _mP, _mP2;
            const double _t_p, _t_m;
            const double _chi_T_pl, _chi_L_pl;
            const double _chi_T_mi, _chi_L_mi;

            static constexpr double chi_1m = 1.0;
            static constexpr double chi_1p = 1.0;
            static constexpr double nI = 2.0;     // [BGL1997] isospin Clebsch-Gordan factor for B->D* and B->D

            double inline _m_c_hat() const { return _m_c_pole() / _m_b_pole(); }

// TODO same function used in <BToDstar>
            double _z(const double & t, const double & t_0) const;
 
// TODO same function used in <BToDstar>
            double _chi_T(const double & u) const;

// TODO same function used in <BToDstar>
            double _chi_L(const double & u) const;

// TODO same function used in <BToDstar>
            double _phi(const double & s, const unsigned & K, const unsigned & a, const unsigned & b, const unsigned & c, const double & chi) const;

            static std::string _par_name(const std::string & ff_name);

        public:

            BGL1997FormFactors(const Parameters &, const Options &);

            ~BGL1997FormFactors();

            static FormFactors<PToP> * make(const Parameters & parameters, const Options & options);

            virtual double f_p(const double & s) const;

            virtual double f_0(const double & s) const;

            virtual double f_t(const double & s) const;
    };
}

#endif