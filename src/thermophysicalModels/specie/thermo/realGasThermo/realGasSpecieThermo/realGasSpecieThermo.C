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

\*---------------------------------------------------------------------------*/

#include "realGasSpecieThermo.H"
#include "IOstreams.H"

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

template<class thermo>
const Foam::scalar Foam::realGasSpecieThermo<thermo>::tol_ = 1.0e-9;

template<class thermo>
const int Foam::realGasSpecieThermo<thermo>::maxIter_ = 500;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
Foam::realGasSpecieThermo<thermo>::realGasSpecieThermo(Istream& is)
:
    thermo(is)
{
    is.check("realGasSpecieThermo::realGasSpecieThermo(Istream& is)");
    
	 
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

template<class thermo>
Foam::Ostream& Foam::operator<<(Ostream& os, const realGasSpecieThermo<thermo>& st)
{
    os  << static_cast<const thermo&>(st);

    os.check("Ostream& operator<<(Ostream& os, const realGasSpecieThermo& st)");
    return os;
}


// ************************************************************************* //
