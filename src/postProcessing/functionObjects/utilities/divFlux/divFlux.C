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

#include "divFlux.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvc.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(divFlux, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        divFlux,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::divFlux::divFlux
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    phiName_(dict.lookup("phiName"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating divFlux for field "
        << phiName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::divFlux::start()
{
    return true;
}


bool Foam::divFlux::execute()
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    if (mesh.foundObject<surfaceScalarField>(phiName_))
    {
        const surfaceScalarField& phi =
            mesh.lookupObject<surfaceScalarField>(phiName_);

        volScalarField divFlux
        (
            IOobject
            (
                "divFlux",
                phi.instance(),
                mesh
            ),
            mag(fvc::div(phi))
        );

        Info<< "Flux divergence min = " << Foam::min(divFlux.internalField())
            << " max = " << Foam::max(divFlux.internalField())
            << " average: " << divFlux.weightedAverage(mesh.V()).value()
            << endl;

        if (mesh.time().outputTime())
        {
            Info << "Writing divFlux field" << endl;
            divFlux.write();
        }

        return true;
    }
    else
    {
        Info<< "Field "  << phiName_ << " not found.  Skipping."
            << endl;

        return false;
    }
}


bool Foam::divFlux::read(const dictionary& dict)
{
    phiName_ = word(dict.lookup("phiName"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}

// ************************************************************************* //
