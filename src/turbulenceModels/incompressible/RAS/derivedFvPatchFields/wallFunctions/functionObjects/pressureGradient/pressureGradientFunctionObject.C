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

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "pressureGradientFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureGradientFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        pressureGradientFunctionObject,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureGradientFunctionObject::pressureGradientFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    regionName_(dict.lookupOrDefault<word>("region", polyMesh::defaultRegion)),
    pName_(dict.lookup("pName")),
    gradp_
    (
        IOobject
        (
            "pressureGradient",
            time_.timeName(),
            time_.lookupObject<fvMesh>(regionName_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        time_.lookupObject<fvMesh>(regionName_),
        // Note: dimensions reset in the constructor body
        dimensionedVector("zero", dimless, vector::zero)
    )
{
    Info<< "Creating pressureGradientFunctionObject for "
        << pName_ << " field." << endl;

    // Reset dimensions
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if (mesh.foundObject<volScalarField>(pName_))
    {
        // Found pressure, reset dimensions for pressure gradient
        const volScalarField& p =
            mesh.lookupObject<volScalarField>(pName_);

        gradp_.dimensions().reset(p.dimensions()/dimLength);
    }
    else
    {
        WarningIn
        (
            "pressureGradientFunctionObject::"
            "pressureGradientFunctionObject"
            "\n("
            "\n    const word& name,"
            "\n    const Time& t,"
            "\n    const dictionary& dict"
            "\n)"
        )   << "Could not find " << pName_ << " field."
            << " Possible dimensions mismatch."
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::pressureGradientFunctionObject::start()
{
    return true;
}


bool Foam::pressureGradientFunctionObject::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    if (mesh.foundObject<volScalarField>(pName_))
    {
        const volScalarField& p =
            mesh.lookupObject<volScalarField>(pName_);

        gradp_ = fvc::grad(p);

        return true;
    }
    else
    {
        InfoIn("bool pressureGradientFunctionObject::execute()")
            << "Field " << pName_ << " not found.  Returning."
            << endl;

        return false;
    }
}


bool Foam::pressureGradientFunctionObject::read(const dictionary& dict)
{
    pName_ = word(dict.lookup("pName"));

    return false;
}


// ************************************************************************* //
