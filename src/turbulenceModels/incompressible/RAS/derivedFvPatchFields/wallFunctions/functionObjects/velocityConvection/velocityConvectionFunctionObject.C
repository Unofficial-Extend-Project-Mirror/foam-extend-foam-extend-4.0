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


bool Foam::velocityConvectionFunctionObject::execute(const bool forceWrite)
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
