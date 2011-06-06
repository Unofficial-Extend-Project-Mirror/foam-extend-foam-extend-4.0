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
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

\*---------------------------------------------------------------------------*/

#include "mixingPlaneCheckFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "mixingPlaneFvsPatchFields.H"
#include "surfaceFields.H"
#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(mixingPlaneCheckFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        mixingPlaneCheckFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mixingPlaneCheckFunctionObject::mixingPlaneCheckFunctionObject
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

    Info << "Creating mixingPlane check" << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::mixingPlaneCheckFunctionObject::start()
{
    return true;
}


bool Foam::mixingPlaneCheckFunctionObject::execute()
{
    const polyMesh& mesh =
        time_.lookupObject<polyMesh>(regionName_);

//     const surfaceScalarField& phi =
//         mesh.lookupObject<surfaceScalarField>(phiName_);

    boolList visited(mesh.boundaryMesh().size(), false);

    forAll (mesh.boundaryMesh(), patchI)
    {
        if (isA<mixingPlanePolyPatch>(mesh.boundaryMesh()[patchI]))
        {
            if (!visited[patchI])
            {
                visited[patchI] = true;

                const mixingPlanePolyPatch& mixingMaster =
                    refCast<const mixingPlanePolyPatch>
                    (
                        mesh.boundaryMesh()[patchI]
                    );

                const label shadowPatchI = mixingMaster.shadowIndex();

                visited[shadowPatchI] = true;

                const mixingPlanePolyPatch& mixingShadow =
                    mixingMaster.shadow();

                // Get access to the mixing plane patch
                const standAlonePatch& mixingPlanePatch =
                    mixingMaster.patchToPatch().mixingPlanePatch();

                // Calculate areas of the mixing patch
                scalarField mixingPlanePatchAreas(mixingPlanePatch.size());
                const vectorField& mixingPlanePatchPoints =
                    mixingPlanePatch.points();

                forAll (mixingPlanePatch, faceI)
                {
                    mixingPlanePatchAreas[faceI] =
                        mixingPlanePatch[faceI].mag(mixingPlanePatchPoints);
                }

                const scalarField masterAreas = mag(mixingMaster.faceAreas());
                const scalarField shadowAreas = mag(mixingShadow.faceAreas());

                scalar sumMasterAreas = gSum(masterAreas);
                scalar sumShadowAreas = gSum(shadowAreas);
                scalar sumMixingAreas = sum(mixingPlanePatchAreas);


                Info<< "Mixing plane area check " << nl
                    << "     "  << "Master " << mixingMaster.name()
                    << " = " << sumMasterAreas << nl
                    << "     "  << "Shadow  " << mixingMaster.shadow().name()
                    << " = " << sumShadowAreas << nl
                    << "     "  << "Mixing = " << sumMixingAreas << nl
                    << "     "  << "mixingPlanePatchAreas = " << mixingPlanePatchAreas
                    //<< "     "  << "mixingPlanePatchPoints = " << mixingPlanePatchPoints
                    << endl;


                // Calculate master to strip sum
                scalarField masterToStripsAreas(mixingPlanePatch.size(), 0);

                const labelListList& mppAddr =
                    mixingMaster.patchToPatch().masterProfileToPatchAddr();

                const scalarListList& mppWeights =
                    mixingMaster.patchToPatch().masterProfileToPatchWeights();

                forAll (masterAreas, masterI)
                {
                    const labelList& curMppAddr = mppAddr[masterI];
                    const scalarList& curMppWeights = mppWeights[masterI];

                    forAll (curMppAddr, i)
                    {
                        masterToStripsAreas[curMppAddr[i]] +=
                            curMppWeights[i]*masterAreas[masterI];
                    }
                }

                Info<< "Master scaling = "
                    << mixingPlanePatchAreas/masterToStripsAreas << endl;

                // Calculate shadow to strip sum
                scalarField shadowToStripsAreas(mixingPlanePatch.size(), 0);

                const labelListList& sppAddr =
                    mixingMaster.patchToPatch().slaveProfileToPatchAddr();

                const scalarListList& sppWeights =
                    mixingMaster.patchToPatch().slaveProfileToPatchWeights();

                forAll (shadowAreas, shadowI)
                {
                    const labelList& curSppAddr = sppAddr[shadowI];
                    const scalarList& curSppWeights = sppWeights[shadowI];

                    forAll (curSppAddr, i)
                    {
                        shadowToStripsAreas[curSppAddr[i]] +=
                            curSppWeights[i]*shadowAreas[shadowI];
                    }
                }

                Info<< "Shadow scaling = "
                    << mixingPlanePatchAreas/shadowToStripsAreas << endl;
            }
        }
    }


//         {
//	    word patchName = phi.boundaryField()[patchI].patch().name();

//             if (patchName == masterPatchName_ && !visited[patchI])
//             {
//                 visited[patchI] = true;

//                 // Calculate local and shadow flux
//                 scalar localFlux    = masterPatchScaleFactor_ * gSum(phi.boundaryField()[patchI]);
//                 //scalar localFluxMag = masterPatchScaleFactor_ * gSumMag(phi.boundaryField()[patchI]);
// 		scalar localFluxMag = mag(localFlux);

//                 const mixingPlanePolyPatch& mixingPlanePatch =
//                     refCast<const mixingPlanePolyPatch>
//                     (
//                         phi.boundaryField()[patchI].patch().patch()
//                     );

//                 const label shadowPatchI = mixingPlanePatch.shadowIndex();

//                 visited[shadowPatchI] = true;

//                 scalar shadowFlux    = shadowPatchScaleFactor_ * gSum(phi.boundaryField()[shadowPatchI]);
//                 //scalar shadowFluxMag = shadowPatchScaleFactor_ * gSumMag(phi.boundaryField()[shadowPatchI]);
// 		scalar shadowFluxMag = mag(shadowFlux);

//                 Info<< "mixingPlane pair " << name_ << " (" << mixingPlanePatch.name() << ", " << mixingPlanePatch.shadow().name() << ") : "
// 		    << " flux: " << localFlux << " " << shadowFlux
// 		    << " : mag: " <<  localFluxMag << " " << shadowFluxMag
//                     << " Diff = " << localFlux + shadowFlux << " or "
//                     << mag(localFlux + shadowFlux)/(localFluxMag + SMALL)*100
//                     << " %" << endl;
//             }
//         }
//     }

    return true;
}


bool Foam::mixingPlaneCheckFunctionObject::read(const dictionary& dict)
{
    return false;
}


// ************************************************************************* //
