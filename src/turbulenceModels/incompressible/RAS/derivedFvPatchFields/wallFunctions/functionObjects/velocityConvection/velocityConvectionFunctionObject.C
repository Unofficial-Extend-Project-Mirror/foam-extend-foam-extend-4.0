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
    Vuko Vukcevic, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "velocityConvectionFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "fvcDiv.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(velocityConvectionFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        velocityConvectionFunctionObject,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::velocityConvectionFunctionObject::velocityConvectionFunctionObject
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
    UName_(dict.lookup("UName")),
    phiName_(dict.lookup("phiName")),
    convection_
    (
        IOobject
        (
            "velocityConvection",
            time_.timeName(),
            time_.lookupObject<fvMesh>(regionName_),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        time_.lookupObject<fvMesh>(regionName_),
        // Note: dimensions will be overwritten on start
        dimensionedVector("zero", dimless, vector::zero)
    )
{
    Info<< "Creating velocityConvectionFunctionObject for "
        << UName_ << " and " << phiName_ << " field." << endl;

    // Reset dimensions
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if
    (
        mesh.foundObject<volVectorField>(UName_)
     && mesh.foundObject<surfaceScalarField>(phiName_)
    )
    {
        // Found velocity and flux, reset dimensions for convection
        const volVectorField& U =
            mesh.lookupObject<volVectorField>(UName_);

        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        convection_.dimensions().reset
        (
            phi.dimensions()*U.dimensions()/dimArea/dimLength
        );
    }
    else
    {
        WarningIn
        (
            "velocityConvectionFunctionObject::"
            "velocityConvectionFunctionObject"
            "\n("
            "\n    const word& name,"
            "\n    const Time& t,"
            "\n    const dictionary& dict"
            "\n)"
        )   << "Could not find " << UName_ << " and " << phiName_
            << " fields. Possible dimensions mismatch."
            << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::velocityConvectionFunctionObject::start()
{
    return true;
}


bool Foam::velocityConvectionFunctionObject::execute()
{
    const fvMesh& mesh = time_.lookupObject<fvMesh>(regionName_);

    if
    (
        mesh.foundObject<volVectorField>(UName_)
     && mesh.foundObject<surfaceScalarField>(phiName_)
    )
    {
        // Found velocity and flux, reset dimensions for convection
        const volVectorField& U =
            mesh.lookupObject<volVectorField>(UName_);

        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        convection_ = fvc::div(phi, U);

        return true;
    }
    else
    {
        InfoIn("bool velocityConvectionFunctionObject::execute()")
            << "Field " << UName_ << " or " << phiName_ << " not found."
            << " Returning."
            << endl;

        return false;
    }
}


bool Foam::velocityConvectionFunctionObject::read(const dictionary& dict)
{
    UName_ = word(dict.lookup("UName"));
    phiName_ = word(dict.lookup("phiName"));

    return false;
}


// ************************************************************************* //
