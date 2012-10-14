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
    Soave Redlich Kwong equation of state for mixtures.

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig 
Germany

\*---------------------------------------------------------------------------*/

#include "mixtureSoaveRedlichKwong.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(Istream& is)
:
    soaveRedlichKwong(is),
    numOfComp(1),
    singleComponent(1)
{ 
    //CL: set size of weigths, mixtureComponents ... to 10,
    //CL: when more mixture componentents are used
    //CL: size of the DynamicLis increases automatically
    weigths.setSize(10);
    mixtureComponents.setSize(10);
    aComponents.setSize(10);
    daComponents.setSize(10);
    d2aComponents.setSize(10);

    //Save a pointer of this object in the mixtureComponents array
    mixtureComponents[0]=this;
    is.check("mixtureSoaveRedlichKwong::mixtureSoaveRedlichKwong(Istream& is)");
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const mixtureSoaveRedlichKwong& srk)
{
    os  << static_cast<const soaveRedlichKwong&>(srk)<< tab;

    os.check("Ostream& operator<<(Ostream& os, const mixtureSoaveRedlichKwong& st)");
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
