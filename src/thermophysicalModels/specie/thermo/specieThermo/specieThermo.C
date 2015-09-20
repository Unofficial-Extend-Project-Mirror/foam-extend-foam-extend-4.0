/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

\*---------------------------------------------------------------------------*/

#include "specieThermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

template<class thermo>
const Foam::debug::tolerancesSwitch
Foam::specieThermo<thermo>::tol_
(
    "speciesThermoTol",
    1e-4
);


template<class thermo>
const Foam::debug::tolerancesSwitch
Foam::specieThermo<thermo>::TJump_
(
    "speciesThermoTJump",
    20
);


template<class thermo>
const Foam::debug::optimisationSwitch
Foam::specieThermo<thermo>::maxIter_
(
    "speciesThermoMaxIter",
    100
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
Foam::specieThermo<thermo>::specieThermo(Istream& is)
:
    thermo(is)
{
    is.check("specieThermo::specieThermo(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const specieThermo<thermo>& st)
{
    os  << static_cast<const thermo&>(st);

    os.check("Ostream& operator<<(Ostream& os, const specieThermo& st)");
    return os;
}


// ************************************************************************* //
