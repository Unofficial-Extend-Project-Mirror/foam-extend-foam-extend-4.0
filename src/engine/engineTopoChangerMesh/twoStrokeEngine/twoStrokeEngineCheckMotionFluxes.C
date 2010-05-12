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

\*---------------------------------------------------------------------------*/


#include "twoStrokeEngine.H"
#include "slidingInterface.H"
#include "layerAdditionRemoval.H"
#include "surfaceFields.H"
#include "regionSplit.H"
#include "componentMixedTetPolyPatchVectorField.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
void Foam::twoStrokeEngine::checkMotionFluxes()
{
/*

    //checking the motion fluxes for the cutFaceZone

    Info << "Checking the motion fluxes in the cut-face zone" << endl;

    const labelList& cutFaceZoneAddressing =
        faceZones()[faceZones().findZoneID("cutFaceZone")];

    boolList calculatedMeshPhi(V().size(), false);
    scalarField sumMeshPhi(V().size(), 0.0);

    forAll(cutFaceZoneAddressing, i)
    {
        label facei = cutFaceZoneAddressing[i];

        calculatedMeshPhi[owner()[facei]] = true;
        calculatedMeshPhi[neighbour()[facei]] = true;
    }

    Info << "Checking the motion fluxes in the liner zone" << endl;

    const labelList& linerAddressing =
        faceZones()[faceZones().findZoneID(scavInCylPatchName_ + "Zone")];

    forAll(linerAddressing, i)
    {
        label facei = linerAddressing[i];

        calculatedMeshPhi[owner()[facei]] = true;
        calculatedMeshPhi[neighbour()[facei]] = true;
    }

    const polyPatch& wallPatch =
        boundaryMesh()[boundaryMesh().findPatchID("wall")];

    labelList wallLabels(wallPatch.size());

    forAll (wallLabels, i)
    {
        wallLabels[i] = wallPatch.start() + i;
    }

    forAll(wallLabels, i)
    {
        label facei = wallLabels[i];

        calculatedMeshPhi[owner()[facei]] = true;
        calculatedMeshPhi[neighbour()[facei]] = true;
    }

    forAll(owner(), facei)
    {
        sumMeshPhi[owner()[facei]] += phi()[facei];
        sumMeshPhi[neighbour()[facei]] -= phi()[facei];
    }

    forAll(boundary(), patchi)
    {
        const unallocLabelList& pFaceCells =
            boundary()[patchi].faceCells();

        const fvsPatchField<scalar>& pssf = phi().boundaryField()[patchi];

        forAll(boundary()[patchi], facei)
        {
            sumMeshPhi[pFaceCells[facei]] += pssf[facei];
        }
    }

    sumMeshPhi /= V();

    DynamicList<label> checkMeshPhi;
    label checkMeshPhiSize = 0;

    forAll(calculatedMeshPhi, i)
    {
        if(calculatedMeshPhi[i] == true)
        {
            checkMeshPhi.append(i);
            checkMeshPhiSize++;
        }
    }

    Info << "checkMeshPhiSize = " << checkMeshPhiSize << endl;

    checkMeshPhi.setSize(checkMeshPhiSize);

    scalarField dVdt(checkMeshPhiSize, 0.0);

    forAll(dVdt, i)
    {
        label celli = checkMeshPhi[i];

        dVdt[i] = (1.0 - V0()[celli]/V()[celli])/engTime().deltaT().value();
    }

    Info << "checking mesh flux in the cutFaces" << endl;

    label nOutCheck = 0;

    forAll(checkMeshPhi, i)
    {
        scalar sumCheck = dVdt[i] - sumMeshPhi[checkMeshPhi[i]];
        if(mag(sumCheck) > 1)
        {
            Info<< "LOCAL sumCheck[" << checkMeshPhi[i] <<  "] = "
                << sumCheck << endl;
            nOutCheck++;
        }
    }

    Info<< "found " << nOutCheck << " cells with inconsistent motion fluxes"
        << endl;
    Info << "end " <<  endl;

*/
    {

        scalarField dVdt = (1.0 - V0()/V())/engTime().deltaT();

        scalarField sumMeshPhi(V().size(), 0.0);

        forAll(owner(), facei)
        {
            sumMeshPhi[owner()[facei]] += phi()[facei];
            sumMeshPhi[neighbour()[facei]] -= phi()[facei];
        }

        forAll(boundary(), patchi)
        {
            const unallocLabelList& pFaceCells =
                boundary()[patchi].faceCells();

            const fvsPatchField<scalar>& pssf = phi().boundaryField()[patchi];

            forAll(boundary()[patchi], facei)
            {
                sumMeshPhi[pFaceCells[facei]] += pssf[facei];
            }
        }

        sumMeshPhi /= V();

        scalarField sumCheck = dVdt - sumMeshPhi;

        label nOutCheck = 0;

        forAll(sumCheck, celli)
        {
            if(mag(sumCheck[celli]) > 1)
            {
                Info<< "sumCheck GLOBAL[" << celli << "] = "
                    << sumCheck[celli] << endl;
                nOutCheck++;
            }
        }

        Info<< "found " << nOutCheck
            << " cells with inconsistent motion fluxes GLOBAL" << endl;

    }
}
