/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Base class of the thermophysical property types.

\*---------------------------------------------------------------------------*/

#include "specie.H"
#include "IOstreams.H"
#include "dimensionedConstants.H"

/* * * * * * * * * * * * * public constants  * * * * * * * * * * * * */

//- Universal gas constant (default in [J/(kmol K)])
const Foam::scalar Foam::specie::RR = dimensionedConstant("R", 8314.51);

//- Standard pressure (default in [Pa])
const Foam::scalar Foam::specie::Pstd = dimensionedConstant("Pstd", 1.0e5);

//- Standard temperature (default in [K])
const Foam::scalar Foam::specie::Tstd = dimensionedConstant("Tstd", 298.15);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::specie::specie(Istream& is)
:
    name_(is),
    nMoles_(readScalar(is)),
    molWeight_(readScalar(is))
{
    is.check("specie::specie(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const specie& st)
{
    os  << st.name_ << tab
        << st.nMoles_ << tab
        << st.molWeight_;

    os.check("Ostream& operator<<(Ostream& os, const specie& st)");
    return os;
}


// ************************************************************************* //
