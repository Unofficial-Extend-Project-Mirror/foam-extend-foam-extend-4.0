/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
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
    a0_(0.457235*pow(this->RR(), 2)*pow(Tcrit_, 2)/pcrit_),
    b_(0.077796*this->RR()*Tcrit_/pcrit_),
    n_(0.37464 + 1.54226*azentricFactor_ - 0.26992*pow(azentricFactor_, 2)),
    b2_(b_*b_),
    //CL: Only uses the default values
    rhoMin_(1e-3),
    rhoMax_(1500),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0),
    // Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(this->Pstd(), this->Tstd(), this->Pstd()/(this->Tstd()*this->R())))
{
    is.check("pengRobinson::pengRobinson(Istream& is)");
}
//CL: Constructed needed in OpenFOAM 2.x.x
//CL: Code works fine, but compiling problem in OpenFOAM 1.6.ext
//CL:  because specie has no constructor using dict
/*
pengRobinson::pengRobinson(const dictionary& dict)
:
    specie(dict),
    pcrit_(readScalar(dict.subDict("equationOfState").lookup("pCritical"))),
    Tcrit_(readScalar(dict.subDict("equationOfState").lookup("TCritical"))),
    azentricFactor_(readScalar(dict.subDict("equationOfState").lookup("azentricFactor"))),
    //CL: rhoMin and rhoMax are only used as boundaries for the bisection methode (see rho function)
    //CL: important: rhoMin and rhoMax are not used as boundary for the newton solver
    //CL: therefore, rho can be larger than rhoMax and smaller than rhoMin
    a0_(0.457235*pow(this->RR(), 2)*pow(Tcrit_, 2)/pcrit_),
    b_(0.077796*this->RR()*Tcrit_/pcrit_),
    n_(0.37464 + 1.54226*azentricFactor_ - 0.26992*pow(azentricFactor_, 2)),
    b2_(b_*b_),
    b3_(b2_*b_),
    b4_(b3_*b_),
    b5_(b4_*b_),
    b6_(b5_*b_),
    rhoMin_(dict.subDict("equationOfState").lookupOrDefault("rhoMin",1e-3)),
    rhoMax_(dict.subDict("equationOfState").lookupOrDefault("rhoMax",1500)),
    aSave(0.0),
    daSave(0.0),
    d2aSave(0.0),
    TSave(0.0),
    // Starting GUESS for the density by ideal gas law
    rhostd_(this->rho(this->Pstd(), this->Tstd(), this->Pstd()/(this->Tstd()*this->R())))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pengRobinson::write(Ostream& os) const
{
    specie::write(os);

    dictionary dict("equationOfState");
    dict.add("pCritical", pcrit_);
    dict.add("TCritical", Tcrit_);
    dict.add("azentricFactor", azentricFactor_);
    dict.add("rhoMin", rhoMin_);
    dict.add("rhoMax", rhoMax_);

    os  << indent << dict.dictName() << dict;
}
*/


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const pengRobinson& pr)
{
    os  << static_cast<const specie&>(pr)<< token::SPACE
        << pr.pcrit_ << tab<< pr.Tcrit_<< tab << pr.azentricFactor_;

    os.check("Ostream& operator<<(Ostream& os, const pengRobinson& st)");
    return os;
}


// ************************************************************************* //
