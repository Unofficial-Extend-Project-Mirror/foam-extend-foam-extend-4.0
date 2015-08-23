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
    Soave Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "soaveRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::soaveRedlichKwong::soaveRedlichKwong(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    azentricFactor_(readScalar(is)),
    a0_(0.42747*pow(this->RR,2)*pow(Tcrit_,2)/(pcrit_)),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
    n_(0.48508+1.55171*azentricFactor_-0.15613*pow(azentricFactor_,2)),
    //CL: Only uses the default values
    rhoMin_(1e-3),
    rhoMax_(1500),
    // Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(Pstd,Tstd,Pstd/(Tstd*this->R()))),
    b2_(b_*b_),
    b3_(b2_*b_),
    b5_(b2_*b3_),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0)
{
    is.check("soaveRedlichKwong::soaveRedlichKwong(Istream& is)"); 	
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const soaveRedlichKwong& srk)
{
    os  << static_cast<const specie&>(srk)<< token::SPACE 
        << srk.pcrit_ << tab<< srk.Tcrit_<<tab<<srk.azentricFactor_;

    os.check("Ostream& operator<<(Ostream& os, const soaveRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
