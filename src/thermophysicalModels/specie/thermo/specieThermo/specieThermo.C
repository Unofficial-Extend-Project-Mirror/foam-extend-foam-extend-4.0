/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Description
    Basic thermodynamics type based on the use of fitting functions for
    cp, h, s obtained from the template argument type thermo.  All other
    properties are derived from these primitive functions.

\*---------------------------------------------------------------------------*/

#include "specieThermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * Private static data * * * * * * * * * * * * * */

template<class thermo>
const Foam::scalar Foam::specieThermo<thermo>::tol_
(
    debug::tolerances("speciesThermoTol", 1e-4)
);


template<class thermo>
const Foam::scalar Foam::specieThermo<thermo>::TJump_
(
    debug::tolerances("speciesThermoTJump", 20)
);


template<class thermo>
const int Foam::specieThermo<thermo>::maxIter_
(
    debug::optimisationSwitch("speciesThermoMaxIter", 100)
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
