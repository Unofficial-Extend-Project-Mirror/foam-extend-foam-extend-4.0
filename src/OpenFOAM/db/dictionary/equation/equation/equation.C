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

\*---------------------------------------------------------------------------*/

#include "equation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::equation::typeName = "equation";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equation::equation()
:
    equationName_(word::null),
    rawText_(""),
    lastResult_(word::null, dimless, 0),
    overrideDimensions_(dimless),
    changeDimensions_(false)
{
//    internalScalars_ = new scalarList(0);
    setSize(0);
}


Foam::equation::equation(Istream& is, const word& name)
:
    equationName_(name),
    rawText_(""),
    lastResult_(word::null, dimless, 0),
    overrideDimensions_(dimless),
    changeDimensions_(false)
{
//    internalScalars_ = new scalarList(0);
    operator>>(is, *this);
}


Foam::equation::equation
(
    const Foam::word& equationName,
    const Foam::string& rawText,
    const Foam::dimensionSet& overrideDimensions,
    const bool& changeDimensions
)
:
    equationName_(equationName),
    rawText_(rawText),
    lastResult_(equationName, dimless, 0),
    overrideDimensions_(overrideDimensions),
    changeDimensions_(changeDimensions)
{
//    internalScalars_ = new scalarList(0);
    setSize(0);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equation::~equation()
{}

// * * * * * * * * * * * * * * * * * Operators * * * * * * * * * * * * * * * //

void Foam::equation::operator=(Foam::equation& eqn)
{
    equationName_ = eqn.equationName_;
    rawText_ = eqn.rawText_;
//    this->operations() = eqn.operations();
    overrideDimensions_.reset(eqn.overrideDimensions_);
    changeDimensions_ = eqn.changeDimensions_;
}



// ************************************************************************* //
