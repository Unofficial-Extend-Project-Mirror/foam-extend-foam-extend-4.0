/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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
    Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany


\*---------------------------------------------------------------------------*/

#include "redlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Private static data * * * * * * * * * * * * * */

const scalar redlichKwong::rhoMin_
(
    debug::tolerances("redlichKwongRhoMin", 1e-3)
);

const scalar redlichKwong::rhoMax_
(
    debug::tolerances("redlichKwongRhoMax", 1500)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

redlichKwong::redlichKwong(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    a_(0.42748*pow(this->RR,2)*pow(Tcrit_,2.5)/pcrit_),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
	// Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(Pstd, Tstd, Pstd*this->W()/(Tstd*this->R())))
{ 
    is.check("redlichKwong::redlichKwong(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const redlichKwong& pg)
{
    os  << static_cast<const specie&>(pg)<< tab
        << pg.pcrit_ << tab<< pg.Tcrit_;

    os.check("Ostream& operator<<(Ostream& os, const redlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
