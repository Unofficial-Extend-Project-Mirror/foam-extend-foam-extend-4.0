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

Description
    Tracing function. Useful for computing the balance of the flux across
    a mixingPlane Interface.

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "mixingPlaneFvsPatchFields.H"


namespace Foam
{

// Temporary hack: useful for tracking balance of phi across interface
void
traceMixingPlaneFlux(
    surfaceScalarField& phi,
    scalar& masterPatchScaleFactor,
    scalar& shadowPatchScaleFactor,
    string  traceId
)
{
    boolList visited(phi.boundaryField().size(), false);

    forAll (phi.boundaryField(), patchI)
    {
        if (isA<mixingPlaneFvsPatchScalarField>(phi.boundaryField()[patchI]))
	{
            word patchName = phi.boundaryField()[patchI].patch().name();

            const mixingPlanePolyPatch& mixingPlanePatch =
                refCast<const mixingPlanePolyPatch>
                (
                    phi.boundaryField()[patchI].patch().patch()
                );

            if (mixingPlanePatch.master() && !visited[patchI])
	    {
                visited[patchI] = true;

                // Calculate local and shadow flux
                scalar localFlux    = masterPatchScaleFactor * gSum(phi.boundaryField()[patchI]);
                //scalar localFluxMag = masterPatchScaleFactor_ * gSumMag(phi.boundaryField()[patchI]);
                scalar localFluxMag = mag(localFlux);


                const label shadowPatchI = mixingPlanePatch.shadowIndex();

                visited[shadowPatchI] = true;

                scalar shadowFlux    = shadowPatchScaleFactor * gSum(phi.boundaryField()[shadowPatchI]);
                //scalar shadowFluxMag = shadowPatchScaleFactor_ * gSumMag(phi.boundaryField()[shadowPatchI]);
                scalar shadowFluxMag = mag(shadowFlux);

                Info<< "====> traceMixingPlaneFlux::: ID: " << traceId << ":  mixingPlane pair (" << mixingPlanePatch.name() << ", " << mixingPlanePatch.shadow().name() << ") : "
                    << " flux: " << localFlux << " " << shadowFlux
                    << " : mag: " <<  localFluxMag << " " << shadowFluxMag
                    << " Diff = " << localFlux + shadowFlux << " or "
                    << mag(localFlux + shadowFlux)/(localFluxMag + SMALL)*100
                    << " %" << endl;
	    }
	}
    }

}

}
