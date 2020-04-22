/* vim: set sw=4 sts=4 et foldmethod=syntax : */

/*
 * Copyright (c) 2017-2020 Danny van Dyk
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

#include <eos/rare-b-decays/nonlocal-correlator.hh>
#include <eos/utils/options-impl.hh>
#include <eos/utils/private_implementation_pattern-impl.hh>

namespace eos
{
    // B -> P
    NonlocalCorrelator<nc::PToP>::~NonlocalCorrelator()
    {
    }

    complex<double>
    NonlocalCorrelator<nc::PToP>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalCorrelator<nc::PToP>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalCorrelator<nc::PToP>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    // B -> V
    NonlocalCorrelator<nc::PToV>::~NonlocalCorrelator()
    {
    }

    complex<double>
    NonlocalCorrelator<nc::PToV>::jpsi_residues_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the J/psi residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalCorrelator<nc::PToV>::psi2s_residues_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the psi(2S) residues has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }

    complex<double>
    NonlocalCorrelator<nc::PToV>::moments_not_implemented() const
    {
        throw InternalError("A NonlocalCorrelator without implementation of the LCSR moments has been erroneously used.");

        // Quiet down compiler warnings.
        return complex<double>(0.0);
    }
}
