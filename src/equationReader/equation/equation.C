/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#include "equation.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::equation::typeName = "equation";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::equation::equation()
:
    equationName_(word::null),
    ops_(0),
    rawText_(""),
    lastResult_(word::null, dimless, 0),
    overrideDimensions_(dimless),
    changeDimensions_(false)
{}


Foam::equation::equation(const equation& newEqn)
:
    equationName_(newEqn.equationName_),
    ops_(0),
    rawText_(newEqn.rawText_),
    lastResult_(word::null, dimless, 0),
    overrideDimensions_(newEqn.overrideDimensions_),
    changeDimensions_(newEqn.changeDimensions_)
{}

Foam::equation::equation(Istream& is, const word& name)
:
    equationName_(name),
    ops_(0),
    rawText_(""),
    lastResult_(word::null, dimless, 0),
    overrideDimensions_(dimless),
    changeDimensions_(false)
{
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
    ops_(0),
    rawText_(rawText),
    lastResult_(equationName, dimless, 0),
    overrideDimensions_(overrideDimensions),
    changeDimensions_(changeDimensions)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::equation::~equation()
{}


// * * * * * * * * * * * * * * Member functions  * * * * * * * * * * * * * * //

void Foam::equation::clear() const
{
    ops_.clear();
}


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
