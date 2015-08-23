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
    Peng Robinson equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "pengRobinson.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pengRobinson::pengRobinson(Istream& is)
:
    specie(is),
    pcrit_(readScalar(is)),
    Tcrit_(readScalar(is)),
    azentricFactor_(readScalar(is)),
    a0_(0.457235*pow(this->RR,2)*pow(Tcrit_,2)/pcrit_),
    b_(0.077796*this->RR*Tcrit_/pcrit_), 
    n_(0.37464+1.54226*azentricFactor_-0.26992*pow(azentricFactor_,2)),
    b2_(b_*b_),
    b3_(b2_*b_),
    b4_(b3_*b_),
    b5_(b4_*b_),
    b6_(b5_*b_),
    //CL: Only uses the default values
    rhoMax_(1500),  
    rhoMin_(1e-3),
    rhostd_(this->rho(Pstd,Tstd,Pstd/(Tstd*this->R()))),	
    aSave(0.0), 
    daSave(0.0), 
    d2aSave(0.0), 
    TSave(0.0)
{
    is.check("pengRobinson::pengRobinson(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pengRobinson& pr)
{
    os  << static_cast<const specie&>(pr)<< token::SPACE 
        << pr.pcrit_ << tab<< pr.Tcrit_<< tab << pr.azentricFactor_;

    os.check("Ostream& operator<<(Ostream& os, const pengRobinson& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
