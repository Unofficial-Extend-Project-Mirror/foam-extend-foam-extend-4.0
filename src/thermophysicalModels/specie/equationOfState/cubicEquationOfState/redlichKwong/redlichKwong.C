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
    Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "redlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::redlichKwong::redlichKwong(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    a_(0.42748*pow(this->RR,2)*pow(Tcrit_,2.5)/pcrit_),
    b_(0.08664*this->RR*Tcrit_/pcrit_),
    //CL: Only uses the default values
    rhoMax_(1500),
    rhoMin_(1e-3),
    b2_(pow(b_,2)),
    b3_(pow(b_,3)),
    b5_(pow(b_,5)),
    // Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(Pstd, Tstd, Pstd/(Tstd*this->R())))
{ 
    is.check("redlichKwong::redlichKwong(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const redlichKwong& rk)
{
    os  << static_cast<const specie&>(rk)<< token::SPACE 
        << rk.pcrit_ << tab<< rk.Tcrit_;

    os.check("Ostream& operator<<(Ostream& os, const redlichKwong& st)");
    return os;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
