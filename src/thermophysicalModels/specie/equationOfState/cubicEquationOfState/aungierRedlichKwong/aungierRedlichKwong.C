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
    Aungier Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "aungierRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::aungierRedlichKwong::aungierRedlichKwong(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    azentricFactor_(readScalar(is)),
    rhocrit_(readScalar(is)),    
    a0_(0.42747*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
    c_(this->RR*Tcrit_/(pcrit_+(a0_/(this->W()/rhocrit_*(this->W()/rhocrit_+b_))))+b_-this->W()/rhocrit_),
    n_(0.4986+1.2735*azentricFactor_+0.4754*pow(azentricFactor_,2)),
    b2_(pow(b_,2)),
    b3_(pow(b_,3)),
    b4_(pow(b_,4)),
    b5_(pow(b_,5)),
    c2_(pow(c_,2)),
    // Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(Pstd,Tstd,Pstd/(Tstd*this->R()))),
    //CL: Only uses the default values
    rhoMax_(1500),
    rhoMin_(1e-3),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0)
{
    is.check("aungierRedlichKwong::aungierRedlichKwong(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const aungierRedlichKwong& ark)
{
    os  << static_cast<const specie&>(ark)<< token::SPACE 
        << ark.pcrit_ << tab<< ark.Tcrit_<< tab<<ark.azentricFactor_<< tab<<ark.rhocrit_;

    os.check("Ostream& operator<<(Ostream& os, const aungierRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// ************************************************************************* //
