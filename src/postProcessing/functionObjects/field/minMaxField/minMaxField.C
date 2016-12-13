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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "minMaxField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(minMaxField, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        minMaxField,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::minMaxField::minMaxField
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    fieldName_(dict.lookup("name"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating minMaxField for field "
        << fieldName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::minMaxField::start()
{
    return true;
}


bool Foam::minMaxField::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    if (mesh.foundObject<volScalarField>(fieldName_))
    {
        const volScalarField& f = mesh.lookupObject<volScalarField>
        (
            fieldName_
        );

        Info<< "Field " << fieldName_ << " min = " << Foam::min(f).value()
            << " max = " << Foam::max(f).value() << endl;

        return true;
    }
    else if (mesh.foundObject<volVectorField>(fieldName_))
    {
        const volVectorField& f = mesh.lookupObject<volVectorField>(fieldName_);

        volScalarField magF = mag(f);

        Info<< "Field " << fieldName_ << " magnitude min = "
            << Foam::min(magF).value()
            << " max = " << Foam::max(magF).value() << endl;

        return true;
    }
    else
    {
        Info<< "Field "  << fieldName_ << " not found.  Skipping." << endl;

        return false;
    }
}


bool Foam::minMaxField::read(const dictionary& dict)
{
    fieldName_ = word(dict.lookup("name"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}

// ************************************************************************* //
