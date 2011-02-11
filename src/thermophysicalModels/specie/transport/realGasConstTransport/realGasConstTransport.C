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

    Modified version of constTransport, modified to be used with realGases

    Constant properties Transport package.  Templated ito a given
    thermodynamics package (needed for thermal conductivity).

Modified by
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "realGasConstTransport.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class thermo>
realGasConstTransport<thermo>::realGasConstTransport(Istream& is)
:
    thermo(is),
    Mu(readScalar(is)),
    rPr(1.0/readScalar(is))
{
    is.check("realGasConstTransport::realGasConstTransport(Istream& is)");
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class thermo>
Ostream& operator<<(Ostream& os, const realGasConstTransport<thermo>& ct)
{
    operator<<(os, static_cast<const thermo&>(ct));
    os << tab << ct.Mu << tab << 1.0/ct.rPr;

    os.check("Ostream& operator<<(Ostream& os, const realGasConstTransport& ct)");

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
