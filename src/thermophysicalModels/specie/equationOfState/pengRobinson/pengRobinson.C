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
    Peng Robinson equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "pengRobinson.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Private static data * * * * * * * * * * * * * */

const scalar pengRobinson::rhoMin_
(
    debug::tolerances("pengRobinsonRhoMin", 1e-3)
);

const scalar pengRobinson::rhoMax_
(
    debug::tolerances("pengRobinsonRhoMax", 1500)
);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

pengRobinson::pengRobinson(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    azentricFactor_(readScalar(is)),
    a_(0.457235*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_),
    b_(0.077796*this->RR*Tcrit_/pcrit_), 
    n_(0.37464+1.54226*azentricFactor_-0.26992*pow(azentricFactor_,2)), 
    rhostd_(this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R())))  	
{
    is.check("pengRobinson::pengRobinson(Istream& is)");
   
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const pengRobinson& pg)
{
    os  << static_cast<const specie&>(pg)<< tab
        << pg.pcrit_ << tab<< pg.Tcrit_<< tab << pg.azentricFactor_;

    os.check("Ostream& operator<<(Ostream& os, const pengRobinson& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
