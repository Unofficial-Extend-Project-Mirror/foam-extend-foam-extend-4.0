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
    Mixture Aungier Redlich Kwong equation of state.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "mixtureAungierRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixtureAungierRedlichKwong::mixtureAungierRedlichKwong(Istream& is)
:
    aungierRedlichKwong(is),
    Vcrit_(this->W()/rhocrit_),
    Zcrit_(pcrit_/((this->R()*Tcrit_*rhocrit_)))
{
    is.check("mixtureAungierRedlichKwong::mixtureAungierRedlichKwong(Istream& is)");
}
//CL: Constructed needed in OpenFOAM 2.x.x
//CL: Code works fine, but compiling problem in OpenFOAM 1.6.ext
//CL:  because specie has no constructor using dict
/*
mixtureAungierRedlichKwong::mixtureAungierRedlichKwong(const dictionary& dict)
:
    aungierRedlichKwong(dict),
    Vcrit_(this->W()/rhocrit_),
    Zcrit_(pcrit_/((this->R()*Tcrit_*rhocrit_)))
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::mixtureAungierRedlichKwong::write(Ostream& os) const
{
    aungierRedlichKwong::write(os);
}
*/
// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const mixtureAungierRedlichKwong& ark)
{
    os  << static_cast<const specie&>(ark)<< token::SPACE;

    os.check("Ostream& operator<<(Ostream& os, const mixtureAungierRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
