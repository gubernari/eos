/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2015, 2016, 2017 Danny van Dyk
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

#include <eos/rare-b-decays/b-to-kstar-ll-base.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-abbbsw2008.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-bcvdv2016.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-bfs2004.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-cffmpsv2015.hh>
#include <eos/rare-b-decays/b-to-kstar-ll-gp2004.hh>
#include <eos/utils/integrate.hh>
#include <eos/utils/integrate-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>
#include <eos/utils/save.hh>

namespace eos
{
    struct BToKstarDilepton::AngularCoefficients
    {
        double j1s, j1c;
        double j2s, j2c;
        double j3;
        double j4;
        double j5;
        double j6s, j6c;
        double j7;
        double j8;
        double j9;

        AngularCoefficients()
        {
        }

        AngularCoefficients(const std::array<double, 12> & a) :
            j1s(a[0]),
            j1c(a[1]),
            j2s(a[2]),
            j2c(a[3]),
            j3(a[4]),
            j4(a[5]),
            j5(a[6]),
            j6s(a[7]),
            j6c(a[8]),
            j7(a[9]),
            j8(a[10]),
            j9(a[11])
        {
        }
    };

    /*!
     * Implementation for the decay @f$\bar{B} \to \bar{K}^* \ell^+ \ell^-@f$.
     */
    template <>
    struct Implementation<BToKstarDilepton>
    {
        std::shared_ptr<BToKstarDilepton::AmplitudeGenerator> amplitude_generator;

        std::shared_ptr<Model> model;

        UsedParameter hbar;
        UsedParameter m_l;
        UsedParameter tau;
        UsedParameter mu;

        Implementation(const Parameters & p, const Options & o, ParameterUser & u) :
            model(Model::make(o.get("model", "WilsonScan"), p, o)),
            hbar(p["hbar"], u),
            m_l(p["mass::" + o.get("l", "mu")], u),
            tau(p["life_time::B_" + o.get("q", "d")], u),
            mu(p["mu"], u)
        {
            std::string tag = o.get("tag", "");
#if 0
            if ("ABBBSW2008" == tag)
            {
                amplitude_generator.reset(new BToKstarDileptonAmplitudes<tag::ABBBSW2008>(p, o));
            }
            else
#endif
            if ("BCvDV2016" == tag)
            {
                amplitude_generator.reset(new BToKstarDileptonAmplitudes<tag::BCvDV2016>(p, o));
            }
            else if ("BFS2004" == tag)
            {
                amplitude_generator.reset(new BToKstarDileptonAmplitudes<tag::BFS2004>(p, o));
            }
            else if ("GP2004" == tag)
            {
                amplitude_generator.reset(new BToKstarDileptonAmplitudes<tag::GP2004>(p, o));
            }
            else
            {
                throw InternalError("BToKstarDilepton: Unknown tag or no valid tag specified (tag = '" + tag + "')!");
            }

            u.uses(*amplitude_generator);
        }

        ~Implementation()
        {
        }

        inline std::array<double, 12> angular_coefficients_array(const BToKstarDilepton::Amplitudes & A, const double & s) const
        {
            // cf. [BHvD2010], p. 26, eqs. (A1)-(A11)
            // cf. [BHvD2012], app B, eqs. (B1)-(B12)
            std::array<double, 12> result;

            double z = 4.0 * power_of<2>(m_l()) / s;
            double y = m_l / std::sqrt(s);
            double beta2 = 1.0 - z;
            double beta = std::sqrt(beta2);

            // j1s
            result[0] = 3.0 / 4.0 * (
                  (2.0 + beta2) / 4.0 * (norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_para_left) + norm(A.a_para_right))
                  + z * real(A.a_perp_left * conj(A.a_perp_right) + A.a_para_left * conj(A.a_para_right))
                  + 4.0 * beta2 * (norm(A.a_long_perp) + norm(A.a_long_para))
                  + 4.0 * (4.0 - 3.0 * beta2) * (norm(A.a_time_perp) + norm(A.a_time_para))
                  + 8.0 * std::sqrt(2.0) * y * real(
                       (A.a_para_left + A.a_para_right)   * conj(A.a_time_para)
                     + (A.a_perp_left + A.a_perp_right) * conj(A.a_time_perp)
                  )
               );
            // j1c
            result[1] = 3.0 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  + z * (norm(A.a_time) + 2.0 * real(A.a_long_left * conj(A.a_long_right)))
                  + beta2 * norm(A.a_scal)
                  + 8.0 * (2.0 - beta2) * norm(A.a_time_long)
                  + 8.0 * beta2 * norm(A.a_para_perp)
                  + 16.0 * y * real((A.a_long_left + A.a_long_right) * conj(A.a_time_long))
               );
            // j2s
            result[2] = 3.0 * beta2 / 16.0 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) + norm(A.a_para_left) + norm(A.a_para_right)
                  - 16.0 * (norm(A.a_time_perp) + norm(A.a_time_para) + norm(A.a_long_perp) + norm(A.a_long_para))
               );
            // j2c
            result[3] = -3.0 * beta2 / 4.0 * (
                  norm(A.a_long_left) + norm(A.a_long_right)
                  - 8.0 * (norm(A.a_time_long) + norm(A.a_para_perp))
               );
            // j3
            result[4] = 3.0 / 8.0 * beta2 * (
                  norm(A.a_perp_left) + norm(A.a_perp_right) - norm(A.a_para_left) - norm(A.a_para_right)
                  + 16.0 * (norm(A.a_time_para) - norm(A.a_time_perp) + norm(A.a_long_para) - norm(A.a_long_perp))
               );
            // j4
            result[5] = 3.0 / (4.0 * std::sqrt(2.0)) * beta2 * real(
                  A.a_long_left * conj(A.a_para_left) + A.a_long_right * conj(A.a_para_right)
                  - 8.0 * std::sqrt(2.0) * (A.a_time_long * conj(A.a_time_para) + A.a_para_perp * conj(A.a_long_para))
               );
            // j5
            result[6] = 3.0 * std::sqrt(2.0) / 4.0 * beta * real(
                  A.a_long_left * conj(A.a_perp_left) - A.a_long_right * conj(A.a_perp_right)
                  - 2.0 * std::sqrt(2.0) * A.a_time_para * conj(A.a_scal)
                  - y * (
                     (A.a_para_left + A.a_para_right) * conj(A.a_scal)
                     + 4.0 * std::sqrt(2.0) * A.a_long_para * conj(A.a_time)
                     - 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_time_perp)
                     - 4.0 * (A.a_perp_left - A.a_perp_right) * conj(A.a_time_long)
                  )
               );
            // j6s
            result[7] = 3.0 / 2.0 * beta * real(
                  A.a_para_left * conj(A.a_perp_left) - A.a_para_right * conj(A.a_perp_right)
                  + 4.0 * std::sqrt(2.0) * y * (
                     (A.a_perp_left - A.a_perp_right) * conj(A.a_time_para)
                     + (A.a_para_left - A.a_para_right) * conj(A.a_time_perp)
                  )
               );
            // j6c
            result[8] = 3.0 * beta * real(
                  2.0 * A.a_time_long * conj(A.a_scal)
                  + y * (
                     (A.a_long_left + A.a_long_right) * conj(A.a_scal)
                     + 4.0 * A.a_para_perp * conj(A.a_time)
                  )
               );
            // j7
            result[9] = 3.0 * std::sqrt(2.0) / 4.0 * beta * imag(
                  A.a_long_left * conj(A.a_para_left) - A.a_long_right * conj(A.a_para_right)
                  + 2.0 * std::sqrt(2.0) * A.a_time_perp * conj(A.a_scal)
                  + y * (
                     (A.a_perp_left + A.a_perp_right) * conj(A.a_scal)
                     + 4.0 * std::sqrt(2.0) * A.a_long_perp * conj(A.a_time)
                     + 4.0 * std::sqrt(2.0) * (A.a_long_left - A.a_long_right) * conj(A.a_time_para)
                     - 4.0 * (A.a_para_left - A.a_para_right) * conj(A.a_time_long)
                  )
               );
            // j8
            result[10] = 3.0 / 4.0 / std::sqrt(2.0) * beta2 * imag(
                  A.a_long_left * conj(A.a_perp_left) + A.a_long_right * conj(A.a_perp_right)
               );
            // j9
            result[11] = 3.0 / 4.0 * beta2 * imag(
                  conj(A.a_para_left) * A.a_perp_left + conj(A.a_para_right) * A.a_perp_right
               );

            return result;
        }

        inline std::array<double, 12> differential_angular_coefficients_array(const double & s) const
        {
            return angular_coefficients_array(amplitude_generator->amplitudes(s), s);
        }

        inline BToKstarDilepton::AngularCoefficients differential_angular_coefficients(const double & s) const
        {
            return BToKstarDilepton::AngularCoefficients(differential_angular_coefficients_array(s));
        }

        BToKstarDilepton::AngularCoefficients integrated_angular_coefficients(const double & s_min, const double & s_max) const
        {
            std::function<std::array<double, 12> (const double &)> integrand =
                    std::bind(&Implementation<BToKstarDilepton>::differential_angular_coefficients_array, this, std::placeholders::_1);
            std::array<double, 12> integrated_angular_coefficients_array = integrate1D(integrand, 64, s_min, s_max);

            return BToKstarDilepton::AngularCoefficients(integrated_angular_coefficients_array);
        }

        inline double decay_width(const BToKstarDilepton::AngularCoefficients & a_c)
        {
            // cf. [BHvD2010], p. 6, eq. (2.7)
            return 2.0 * a_c.j1s + a_c.j1c - 1.0 / 3.0 * (2.0 * a_c.j2s + a_c.j2c);
        }

        inline double beta_l(const double & s) const
        {
            return std::sqrt(1.0 - 4.0 * m_l * m_l / s);
        }

        double a_fb_zero_crossing() const
        {
            // We trust QCDF results in a validity range from 0.5 GeV^2 < s < 6.0 GeV^2
            static const double min_result = 0.5;
            static const double max_result = 7.0;

            // Use calT_perp / xi_perp = C_7 as start point.
            // Use hard coded values for mu=4.2 GeV and M_B = 5.2795 GeV here.
            WilsonCoefficients<BToS> wc = model->wilson_coefficients_b_to_s(mu(), amplitude_generator->lepton_flavour, amplitude_generator->cp_conjugate);
            const double start = -2.0 * model->m_b_msbar(4.2) * 5.2795 * real(wc.c7() / wc.c9());

            double result = start;
            // clamp result to QCDF validity region
            result = std::max(min_result, result);
            result = std::min(max_result, result);

            // perform a couple of Newton-Raphson steps
            for (unsigned i = 0 ; i < 100 ; ++i)
            {
                double xplus = result * 1.03;
                double xminus = result * 0.97;

                auto a_c_central = differential_angular_coefficients(result);
                double f = a_c_central.j6s + 0.5 * a_c_central.j6c;
                auto a_c_minus   = differential_angular_coefficients(xminus);
                double f_xminus = a_c_minus.j6s + 0.5 * a_c_minus.j6c;
                auto a_c_plus    = differential_angular_coefficients(xplus);
                double f_xplus = a_c_plus.j6s + 0.5 * a_c_plus.j6c;

                double fprime = (f_xplus - f_xminus) / (xplus - xminus);

                if (std::abs(f / fprime) < 1e-8)
                    break;

                result = result - f / fprime;
                // clamp result to QCDF validity region
                result = std::max(min_result, result);
                result = std::min(max_result, result);
            }

            return result;
        }
    };

    BToKstarDilepton::BToKstarDilepton(const Parameters & parameters, const Options & options) :
        PrivateImplementationPattern<BToKstarDilepton>(new Implementation<BToKstarDilepton>(parameters, options, *this))
    {
    }

    BToKstarDilepton::~BToKstarDilepton()
    {
    }

    double
    BToKstarDilepton::a_fb_zero_crossing() const
    {
        return _imp->a_fb_zero_crossing();
    }

    double
    BToKstarDilepton::decay_width(const double & s, const double & c_theta_l, const double & c_theta_k, const double & phi) const
    {
        // compute d^4 Gamma, cf. [BHvD2010], p. 5, Eq. (2.6)
        // Cosine squared of the angles
        double c_theta_k_2 = c_theta_k * c_theta_k;
        double c_theta_l_2 = c_theta_l * c_theta_l;
        double c_phi = cos(phi);
        // Sine squared of the angles
        double s_theta_k_2 = 1.0 - c_theta_k_2;
        double s_theta_l_2 = 1.0 - c_theta_l_2;
        // Sine of the angles
        double s_theta_k = sqrt(s_theta_k_2);
        double s_theta_l = sqrt(s_theta_l_2);
        double s_phi = sin(phi);
        // Cosine of twice the angle
        //double c_2_theta_k = 2.0 * c_theta_k_2 - 1.0;
        double c_2_theta_l = 2.0 * c_theta_l_2 - 1.0;
        double c_2_phi = cos(2.0 * phi);
        // Sine of twice the angle
        double s_2_theta_k = 2.0 * s_theta_k * c_theta_k;
        double s_2_theta_l = 2.0 * s_theta_l * c_theta_l;
        double s_2_phi = sin(2.0 * phi);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        double Gamma = _imp->decay_width(_imp->integrated_angular_coefficients(1.00, 6.00));

        double result = 3.0 / 8.0 / M_PI * (
                 a_c.j1s + (a_c.j1c - a_c.j1s) * c_theta_k_2
                +  (a_c.j2s + (a_c.j2c - a_c.j2s) * c_theta_k_2) * c_2_theta_l
                +  a_c.j3 * s_theta_k_2 * s_theta_l_2 * c_2_phi
                +  a_c.j4 * s_2_theta_k * s_2_theta_l * c_phi
                +  a_c.j5 * s_2_theta_k * s_theta_l * c_phi
                +  (a_c.j6s * s_theta_k_2 + a_c.j6c * c_theta_k_2) * c_theta_l
                +  a_c.j7 * s_2_theta_k * s_theta_l * s_phi
                +  a_c.j8 * s_2_theta_k * s_2_theta_l * s_phi
                +  a_c.j9 * s_theta_k_2 * s_theta_l_2 * s_2_phi
                ) / Gamma;

        return result;
    }
    double
    BToKstarDilepton::differential_decay_width(const double & s) const
    {
        return _imp->decay_width(_imp->differential_angular_coefficients(s));
    }

    double
    BToKstarDilepton::differential_branching_ratio(const double & s) const
    {
        return differential_decay_width(s) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton::differential_forward_backward_asymmetry(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.8)
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j6s + 0.5 * a_c.j6c) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::differential_longitudinal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return (a_c.j1c - a_c.j2c / 3.0) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::differential_transversal_polarisation(const double & s) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_2(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.10)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_3(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.11)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_4(const double & s) const
    {
        // cf. [BHvD2010], p. 6, eq. (2.12)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return std::sqrt((power_of<2>(_imp->beta_l(s) * a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(_imp->beta_l(s) * a_c.j7)));
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_re(const double & s) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.25 * _imp->beta_l(s) * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton::differential_transverse_asymmetry_im(const double & s) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton::differential_p_prime_4(const double & s) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton::differential_p_prime_5(const double & s) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton::differential_p_prime_6(const double & s) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->differential_angular_coefficients(s);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton::differential_h_1(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton::differential_h_2(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton::differential_h_3(const double & s) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return _imp->beta_l(s) * a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton::differential_h_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton::differential_h_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // differential angular coefficients
    double
    BToKstarDilepton::differential_j_1c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1c;
    }

    double
    BToKstarDilepton::differential_j_1s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j1s;
    }

    double
    BToKstarDilepton::differential_j_2c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2c;
    }

    double
    BToKstarDilepton::differential_j_2s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j2s;
    }

    double
    BToKstarDilepton::differential_j_3(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j3;
    }

    double
    BToKstarDilepton::differential_j_4(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j4;
    }

    double
    BToKstarDilepton::differential_j_5(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j5;
    }

    double
    BToKstarDilepton::differential_j_6c(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6c;
    }

    double
    BToKstarDilepton::differential_j_6s(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j6s;
    }

    double
    BToKstarDilepton::differential_j_7(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j7;
    }

    double
    BToKstarDilepton::differential_j_8(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j8;
    }

    double
    BToKstarDilepton::differential_j_9(const double & s) const
    {
        AngularCoefficients a_c = _imp->differential_angular_coefficients(s);
        return a_c.j9;
    }

    double
    BToKstarDilepton::integrated_decay_width(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::integrated_branching_ratio(const double & s_min, const double & s_max) const
    {
        return integrated_decay_width(s_min, s_max) * _imp->tau() / _imp->hbar();
    }

    double
    BToKstarDilepton::integrated_forward_backward_asymmetry(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.8), p. 6
        // cf. [BHvD2012], eq. (A7)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j6s + 0.5 * a_c.j6c) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::integrated_longitudinal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A9)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return (a_c.j1c - a_c.j2c / 3.0) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::integrated_transversal_polarisation(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2012], eq. (A10)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 2.0 * (a_c.j1s - a_c.j2s / 3.0) / _imp->decay_width(a_c);
    }

    double
    BToKstarDilepton::integrated_cp_asymmetry(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        auto gamma     = _imp->decay_width(_imp->integrated_angular_coefficients(s_min, s_max));
        _imp->amplitude_generator->cp_conjugate = true;
        auto gamma_bar = _imp->decay_width(_imp->integrated_angular_coefficients(s_min, s_max));

        return (gamma - gamma_bar) / (gamma + gamma_bar);
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.10), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j3 / a_c.j2s;
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.11), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)) / (-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3)));
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_4(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], eq. (2.12), p. 6
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        return sqrt((power_of<2>(a_c.j5) + 4.0 * power_of<2>(a_c.j8)) / (4.0 * power_of<2>(a_c.j4) + power_of<2>(a_c.j7)));
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [BS2011], eq. (34), p. 9 for the massless case
        return std::sqrt(16.0 * power_of<2>(a_c.j2s) - power_of<2>(a_c.j6s) - 4.0 * (power_of<2>(a_c.j3) + power_of<2>(a_c.j9)))
            / 8.0 / a_c.j2s;
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_re(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (38), p. 10
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.25 * a_c.j6s / a_c.j2s;
    }

    double
    BToKstarDilepton::integrated_transverse_asymmetry_im(const double & s_min, const double & s_max) const
    {
        // cf. [BS2011], eq. (30), p. 8
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return 0.5 * a_c.j9 / a_c.j2s;
    }

    double
    BToKstarDilepton::integrated_p_prime_4(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (15)
        return (a_c.j4 + a_c_bar.j4) / std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s));
    }

    double
    BToKstarDilepton::integrated_p_prime_5(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (16)
        return (a_c.j5 + a_c_bar.j5) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton::integrated_p_prime_6(const double & s_min, const double & s_max) const
    {
        Save<bool> save(_imp->amplitude_generator->cp_conjugate, false);

        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        // cf. [DMRV2012], p. 9, eq. (17)
        return -1.0 * (a_c.j7 + a_c_bar.j7) / (2.0 * std::sqrt(-1.0 * (a_c.j2c + a_c_bar.j2c) * (a_c.j2s + a_c_bar.j2s)));
    }

    double
    BToKstarDilepton::integrated_h_1(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.13)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j4 / sqrt(-a_c.j2c * (2.0 * a_c.j2s - a_c.j3));
    }

    double
    BToKstarDilepton::integrated_h_2(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.14)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return  a_c.j5 / sqrt(-2.0 * a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton::integrated_h_3(const double & s_min, const double & s_max) const
    {
        // cf. [BHvD2010], p. 7, eq. (2.15)
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s / (2.0 * sqrt(power_of<2>(2.0 * a_c.j2s) - power_of<2>(a_c.j3)));
    }

    double
    BToKstarDilepton::integrated_h_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return sqrt(2.0) * a_c.j8 / sqrt(-a_c.j2c * (2.0 * a_c.j2s + a_c.j3));
    }

    double
    BToKstarDilepton::integrated_h_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return -a_c.j9 / sqrt(power_of<2>(2.0 * a_c.j2s) + power_of<2>(a_c.j3));
    }

    // integrated angular coefficients
    double
    BToKstarDilepton::integrated_j_1c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1c;
    }

    double
    BToKstarDilepton::integrated_j_1s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j1s;
    }

    double
    BToKstarDilepton::integrated_j_2c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2c;
    }

    double
    BToKstarDilepton::integrated_j_2s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j2s;
    }

    double
    BToKstarDilepton::integrated_j_3(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j3;
    }

    double
    BToKstarDilepton::integrated_j_4(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j4;
    }

    double
    BToKstarDilepton::integrated_j_5(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j5;
    }

    double
    BToKstarDilepton::integrated_j_6c(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6c;
    }

    double
    BToKstarDilepton::integrated_j_6s(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j6s;
    }

    double
    BToKstarDilepton::integrated_j_7(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j7;
    }

    double
    BToKstarDilepton::integrated_j_8(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j8;
    }

    double
    BToKstarDilepton::integrated_j_9(const double & s_min, const double & s_max) const
    {
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        return a_c.j9;
    }

    // integrated CP-symmetrized angular coefficients (S_i)
    double
    BToKstarDilepton::integrated_s_1s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j1s + a_c_bar.j1s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_1c(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j1c + a_c_bar.j1c) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_2s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j2s + a_c_bar.j2s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_2c(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j2c + a_c_bar.j2c) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_3(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j3 + a_c_bar.j3) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_4(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j4 + a_c_bar.j4) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_5(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j5 + a_c_bar.j5) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_6s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j6s + a_c_bar.j6s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_6c(const double & s_min, const double & s_max) const                    
    {                                                                                                      
      Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);                           
      AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);                     
      _imp->amplitude_generator->cp_conjugate = true;                                                    
      AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);                 
                                                                                                       
      return 4.0 / 3.0 * (a_c.j6c + a_c_bar.j6c) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }                                                                                                      


    double
    BToKstarDilepton::integrated_s_7(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j7 + a_c_bar.j7) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_8(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j8 + a_c_bar.j8) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_s_9(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j9 + a_c_bar.j9) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    // integrated CP-antisymmetrized angular coefficients (A_i)
    double
    BToKstarDilepton::integrated_a_1s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j1s - a_c_bar.j1s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_1c(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j1c - a_c_bar.j1c) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_2s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j2s - a_c_bar.j2s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_2c(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j2c - a_c_bar.j2c) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_3(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j3 - a_c_bar.j3) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_4(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j4 - a_c_bar.j4) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_5(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j5 - a_c_bar.j5) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_6s(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j6s - a_c_bar.j6s) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_7(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j7 - a_c_bar.j7) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_8(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j8 - a_c_bar.j8) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::integrated_a_9(const double & s_min, const double & s_max) const
    {
        Save<bool> cp_conjugate(_imp->amplitude_generator->cp_conjugate, false);
        AngularCoefficients a_c = _imp->integrated_angular_coefficients(s_min, s_max);
        _imp->amplitude_generator->cp_conjugate = true;
        AngularCoefficients a_c_bar = _imp->integrated_angular_coefficients(s_min, s_max);

        return 4.0 / 3.0 * (a_c.j9 - a_c_bar.j9) / (_imp->decay_width(a_c) + _imp->decay_width(a_c_bar));
    }

    double
    BToKstarDilepton::differential_symrel_le_a1v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return pow(ag.m_B() + ag.m_Kstar(), 2) / (2.0 * ag.m_B() * ag.energy(q2))
            * ag.form_factors->a_1(q2) / ag.form_factors->v(q2);
    }

    double
    BToKstarDilepton::differential_symrel_le_t1v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return (ag.m_B() + ag.m_Kstar()) / ag.m_B()
            * ag.form_factors->t_1(q2) / ag.form_factors->v(q2);
    }

    double
    BToKstarDilepton::differential_symrel_le_t2v(const double & q2) const
    {
        const auto & ag = *_imp->amplitude_generator;
        return (ag.m_B() + ag.m_Kstar()) / (2.0 * ag.energy(q2))
            * ag.form_factors->t_2(q2) / ag.form_factors->v(q2);
    }

    const std::string
    BToKstarDilepton::description = "\
The decay Bbar->Kbar^*(-> Kbar pi) l^+ l^-, with l=e,mu,tau \
a charged lepton. Various theory models can be selected using \
the 'tag' option";

    const std::string
    BToKstarDilepton::kinematics_description_s = "\
The invariant mass of the charged lepton pair in GeV^2.";

    const std::string
    BToKstarDilepton::kinematics_description_c_theta_l = "\
The cosine of the negatively-charged lepton l^-'s helicity angle theta_l in the l^+l^- rest frame.";

    const std::string
    BToKstarDilepton::kinematics_description_c_theta_k = "\
The cosine of the Kbar's helicity angle theta_k in the Kbar-pi rest frame.";

    const std::string
    BToKstarDilepton::kinematics_description_phi = "\
The azimuthal angle between the Kbar-pi plane and the l^+l^- plane.";

    /*
     * For diagnostic purposes only!
     */
    BToKstarDilepton::Amplitudes
    BToKstarDilepton::amplitudes(const double & q2) const
    {
        return _imp->amplitude_generator->amplitudes(q2);
    }
}
