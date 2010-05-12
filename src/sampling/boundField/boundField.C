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

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "boundField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "bound.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(boundField, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        boundField,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::boundField::boundField
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(),
    name_(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    fieldName_(dict.lookup("name")),
    minValue_(readScalar(dict.lookup("minValue"))),
    maxValue_(readScalar(dict.lookup("maxValue")))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating boundField for field "
        << fieldName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::boundField::start()
{
    return true;
}


bool Foam::boundField::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    volScalarField& f =
        const_cast<volScalarField&>
        (
            mesh.lookupObject<volScalarField>(fieldName_)
        );

    boundMinMax
    (
        f,
        dimensionedScalar("v", f.dimensions(), minValue_),
        dimensionedScalar("V", f.dimensions(), maxValue_)
    );

    return true;
}


bool Foam::boundField::read(const dictionary& dict)
{
    fieldName_ = word(dict.lookup("name"));
    minValue_ = readScalar(dict.lookup("minValue"));
    maxValue_ = readScalar(dict.lookup("maxValue"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}


// ************************************************************************* //
