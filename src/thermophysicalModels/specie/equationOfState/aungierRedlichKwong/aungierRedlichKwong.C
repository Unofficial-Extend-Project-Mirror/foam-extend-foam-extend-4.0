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
    Aungier Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany



\*---------------------------------------------------------------------------*/

#include "aungierRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

aungierRedlichKwong::aungierRedlichKwong(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
	azentricFactor_(readScalar(is)),
	rhocrit_(readScalar(is))
{
    is.check("aungierRedlichKwong::aungierRedlichKwong(Istream& is)");
    rhostd_=this->rho(Pstd,Tstd,Pstd*this->W()/(Tstd*this->R()));
    rhoMax_=1500;	
    rhoMin_=0.001;
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const aungierRedlichKwong& pg)
{
    os  << static_cast<const specie&>(pg)<< tab
        << pg.pcrit_ << tab<< pg.Tcrit_<<tab<<pg.azentricFactor_<<tab<<pg.rhocrit_;
    os.check("Ostream& operator<<(Ostream& os, const aungierRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
