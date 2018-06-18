/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "FUNCTIONOBJECT.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(FUNCTIONOBJECT, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FUNCTIONOBJECT::FUNCTIONOBJECT
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    obr_(obr),
    wordData_(dict.lookupOrDefault<word>("wordData", "defaultWord")),
    scalarData_(readScalar(dict.lookup("scalarData"))),
    labelData_(readLabel(dict.lookup("labelData")))
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::FUNCTIONOBJECT::~FUNCTIONOBJECT()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::FUNCTIONOBJECT::read(const dictionary& dict)
{
    dict.readIfPresent("wordData", wordData_);
    dict.lookup("scalarData") >> scalarData_;
    dict.lookup("labelData") >> labelData_;
}


void Foam::FUNCTIONOBJECT::execute()
{
    // Do nothing - only valid on write
}


void Foam::FUNCTIONOBJECT::end()
{
    // Do nothing - only valid on write
}


void Foam::FUNCTIONOBJECT::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::FUNCTIONOBJECT::write()
{
}


// ************************************************************************* //
