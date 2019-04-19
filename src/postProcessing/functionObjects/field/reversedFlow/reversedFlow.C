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
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "reversedFlow.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(reversedFlow, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        reversedFlow,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::reversedFlow::reversedFlow
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    patchName_(dict.lookup("patch")),
    fluxFieldName_(dict.lookupOrDefault<word>("flux", "phi"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info<< "Creating reversedFlow for patch "
        << patchName_ << " with flux field "<< fluxFieldName_ << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::reversedFlow::start()
{
    return true;
}


bool Foam::reversedFlow::execute(const bool forceWrite)
{
    const fvMesh& mesh =
        time_.lookupObject<fvMesh>(regionName_);

    const label patchIndex = mesh.boundaryMesh().findPatchID(patchName_);

    if
    (
        mesh.foundObject<surfaceScalarField>(fluxFieldName_)
     && patchIndex > -1
    )
    {
        const surfaceScalarField& phi = mesh.lookupObject<surfaceScalarField>
        (
            fluxFieldName_
        );

        const scalarField& patchPhi = phi.boundaryField()[patchIndex];

        scalar minFlux = gMin(patchPhi);
        scalar netFlux = gSum(patchPhi);
        scalar inFlux = gSum(neg(patchPhi)*patchPhi);
        scalar outFlux = gSum(pos(patchPhi)*patchPhi);

        if (minFlux < -SMALL)
        {
            Info<< "Negative patch flux for patch " << patchName_
                << ".  flux: " << netFlux
                << " in = " << inFlux
                << " out = " << outFlux
                << endl;
        }

        return true;
    }
    else
    {
        Info<< "Patch " << patchName_ << " or flux " << fluxFieldName_
            << " not found.  Skipping." << endl;

        return false;
    }
}


bool Foam::reversedFlow::read(const dictionary& dict)
{
    patchName_ = word(dict.lookup("patch"));
    fluxFieldName_ = word(dict.lookupOrDefault<word>("flux", "phi"));

    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    return false;
}


// ************************************************************************* //
