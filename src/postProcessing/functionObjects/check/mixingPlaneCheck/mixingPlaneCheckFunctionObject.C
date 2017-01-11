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

    Info << "Creating mixingPlane check functionObject" << endl;
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

    boolList visited(mesh.boundaryMesh().size(), false);

    forAll (mesh.boundaryMesh(), patchI)
    {
        if
        (
            isA<mixingPlanePolyPatch>(mesh.boundaryMesh()[patchI])
          && mesh.boundaryMesh()[patchI].size()
        )
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

                // Until the mixingPlane is fully parallelized, we stick with
                // the serial version of sum. The interface is residing on a
                // single processor when running in parallel
                scalar sumMasterAreas = sum(masterAreas);
                scalar sumShadowAreas = sum(shadowAreas);
//                 scalar sumMixingAreas = sum(mixingPlanePatchAreas);

#if 0   // Remove this for now
                Info<< "Mixing plane functionObject: area check " << nl
                    << "     "  << "Master " << mixingMaster.name()
                    << " = " << sumMasterAreas << nl
                    << "     "  << "Shadow  " << mixingMaster.shadow().name()
                    << " = " << sumShadowAreas << nl
                    << "     "  << "Mixing = " << sumMixingAreas << nl
                    << "     "  << "mixingPlanePatchAreas = " << mixingPlanePatchAreas
                    //<< "     "  << "mixingPlanePatchPoints = " << mixingPlanePatchPoints
                    << endl;
#endif

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

#if 0
                // Does not work when in cylindrical coordinates
                Info<< "Master scaling = "
                    << mixingPlanePatchAreas/masterToStripsAreas << endl;
#endif

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

#if 0
                // Does not work when in cylindrical coordinates
                Info<< "Shadow scaling = "
                    << mixingPlanePatchAreas/shadowToStripsAreas << endl;
#endif

               // Old way of computing phi balance

                if (mesh.foundObject<surfaceScalarField>(phiName_))
                {
                    const surfaceScalarField& phi =
                        mesh.lookupObject<surfaceScalarField>(phiName_);

                    // Calculate local and shadow flux
                    scalar masterPatchScaleFactor_ = 1.0;
                    scalar shadowPatchScaleFactor_ = sumMasterAreas/sumShadowAreas;
                    scalar localFlux =masterPatchScaleFactor_*
                        sum(phi.boundaryField()[patchI]);

                    scalar localFluxMag = mag(localFlux);

                    scalar shadowFlux =shadowPatchScaleFactor_*
                        sum(phi.boundaryField()[shadowPatchI]);
//                     scalar shadowFluxMag = mag(shadowFlux);

                    Info<< "Mixing plane pair "
                        << "(" << mixingMaster.name() << ", "
                        << mixingMaster.shadow().name() << ") : "
                        << localFlux << " " << shadowFlux
                        << " Diff = " << localFlux + shadowFlux << " or "
                        << mag(localFlux + shadowFlux)/(localFluxMag + SMALL)*100
                        << " %" << endl;
                }
            }
        }
    }

    return true;
}


bool Foam::mixingPlaneCheckFunctionObject::read(const dictionary& dict)
{
    return false;
}


// ************************************************************************* //
