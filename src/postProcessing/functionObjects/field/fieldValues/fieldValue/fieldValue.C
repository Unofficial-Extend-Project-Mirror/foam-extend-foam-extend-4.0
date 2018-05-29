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

#include "fieldValue.H"
#include "fvMesh.H"
#include "foamTime.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(fieldValue, 0);
    defineRunTimeSelectionTable(fieldValue, dictionary);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::fieldValue::read(const dictionary& dict)
{
    if (active_)
    {
        dict_ = dict;

        log_ = dict.lookupOrDefault<Switch>("log", true);
        dict.lookup("fields") >> fields_;
        dict.lookup("valueOutput") >> valueOutput_;
    }
}


void Foam::fieldValue::write()
{
    if (active_)
    {
        functionObjectFile::write();

        if (log_) Info<< type() << " " << name_ << " output:" << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldValue::fieldValue
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const word& valueType,
    const bool loadFromFiles
)
:
    functionObjectFile(obr, name, valueType),
    name_(name),
    obr_(obr),
    dict_(dict),
    active_(true),
    log_(true),
    sourceName_(word::null),
    fields_(dict.lookup("fields")),
    valueOutput_(dict.lookup("valueOutput")),
    resultDict_(fileName("name"), dictionary::null)
{
    // Only active if obr is an fvMesh
    if (isA<fvMesh>(obr_))
    {
        read(dict);
    }
    else
    {
        WarningIn
        (
            "fieldValue::fieldValue"
            "("
                "const word&, "
                "const objectRegistry&, "
                "const dictionary&, "
                "const bool"
            ")"
        )   << "No fvMesh available, deactivating " << name << nl
            << endl;
        active_ = false;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldValue::~fieldValue()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::fieldValue::execute()
{
    // Do nothing
}


void Foam::fieldValue::end()
{
    // Do nothing
}


void Foam::fieldValue::timeSet()
{
    // Do nothing
}


void Foam::fieldValue::updateMesh(const mapPolyMesh&)
{
    // Do nothing
}


void Foam::fieldValue::movePoints(const pointField&)
{
    // Do nothing
}


// ************************************************************************* //
