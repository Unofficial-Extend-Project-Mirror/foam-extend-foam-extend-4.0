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

\*---------------------------------------------------------------------------*/

#include "ggiCheckFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "ggiFvsPatchFields.H"
#include "cyclicGgiFvsPatchFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(ggiCheckFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        ggiCheckFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ggiCheckFunctionObject::ggiCheckFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    time_(t),
    regionName_(polyMesh::defaultRegion),
    phiName_(dict.lookup("phi"))
{
    if (dict.found("region"))
    {
        dict.lookup("region") >> regionName_;
    }

    Info << "Creating ggi check" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::ggiCheckFunctionObject::start()
{
    return true;
}


bool Foam::ggiCheckFunctionObject::execute()
{
    const objectRegistry& mesh =
        time_.lookupObject<objectRegistry>(regionName_);

    const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

    boolList visited(phi.boundaryField().size(), false);

    forAll (phi.boundaryField(), patchI)
    {
        if (isA<ggiFvsPatchScalarField>(phi.boundaryField()[patchI]))
        {
            if (!visited[patchI])
            {
                visited[patchI] = true;

                // Calculate local and shadow flux
                scalar localFlux = gSum(phi.boundaryField()[patchI]);
                scalar localFluxMag = gSumMag(phi.boundaryField()[patchI]);

                const ggiPolyPatch& ggiPatch =
                    refCast<const ggiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                const label shadowPatchI = ggiPatch.shadowIndex();

                visited[shadowPatchI] = true;

                scalar shadowFlux = gSum(phi.boundaryField()[shadowPatchI]);
                scalar shadowFluxMag =
                    gSumMag(phi.boundaryField()[shadowPatchI]);

                Info<< "GGI pair (" << ggiPatch.name() << ", "
                    << ggiPatch.shadow().name()
                    << ") : " << localFluxMag << " " << shadowFluxMag
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/(localFluxMag + SMALL)*100
                    << " %" << endl;
            }
        }
        else if (isA<cyclicGgiFvsPatchScalarField>(phi.boundaryField()[patchI]))
        {
            if (!visited[patchI])
            {
                visited[patchI] = true;

                // Calculate local and shadow flux
                scalar localFlux = gSum(phi.boundaryField()[patchI]);
                scalar localFluxMag = gSumMag(phi.boundaryField()[patchI]);

                const cyclicGgiPolyPatch& ggiPatch =
                    refCast<const cyclicGgiPolyPatch>
                    (
                        phi.boundaryField()[patchI].patch().patch()
                    );

                const label shadowPatchI = ggiPatch.shadowIndex();

                visited[shadowPatchI] = true;

                scalar shadowFlux = gSum(phi.boundaryField()[shadowPatchI]);
                scalar shadowFluxMag =
                    gSumMag(phi.boundaryField()[shadowPatchI]);

                Info<< "Cyclic GGI pair (" << ggiPatch.name() << ", "
                    << ggiPatch.shadow().name()
                    << ") : " << localFluxMag << " " << shadowFluxMag
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/(localFluxMag + SMALL)*100
                    << " %" << endl;
            }
        }
    }

    return true;
}


bool Foam::ggiCheckFunctionObject::read(const dictionary& dict)
{
    return false;
}

// ************************************************************************* //
