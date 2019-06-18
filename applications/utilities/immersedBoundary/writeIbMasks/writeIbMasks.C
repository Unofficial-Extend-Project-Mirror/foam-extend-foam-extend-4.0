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

Application
    writeIbMasks

Description
    Calculate and write immersed boundary masks

\*---------------------------------------------------------------------------*/

#include "calc.H"
#include "fvc.H"
#include "fvMatrices.H"
#include "immersedBoundaryFvPatch.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::calc(const argList& args, const Time& runTime, const fvMesh& mesh)
{
    Info<< nl << "Calculating gamma" << endl;
    volScalarField gamma
    (
        IOobject
        (
            "gamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 1)
    );
    gamma.internalField() = mesh.V()/mesh.cellVolumes();

    // Report minimal live cell volume
    scalar minLiveGamma = GREAT;
    label minLiveCell = -1;
    const scalarField& gammaIn = gamma.internalField();

    // Collect dead cells
    labelHashSet deadCellsHash;
    
    forAll (mesh.boundary(), patchI)
    {
        if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            const immersedBoundaryFvPatch& ibPatch =
                refCast<const immersedBoundaryFvPatch>
                (
                    mesh.boundary()[patchI]
                );

            const labelList& ibCells = ibPatch.ibPolyPatch().ibCells();

            forAll (ibCells, dcI)
            {
                if (gammaIn[ibCells[dcI]] < minLiveGamma)
                {
                    minLiveGamma = gammaIn[ibCells[dcI]];
                    minLiveCell = ibCells[dcI];
                }
            }

            // Collect dead cells
            deadCellsHash.insert(ibPatch.ibPolyPatch().deadCells());
        }
    }

    Info<< "Min live cell " << minLiveCell
        << " gamma = " << minLiveGamma
        << endl;

    Info<< nl << "Calculating sGamma" << endl;
    surfaceScalarField sGamma
    (
        IOobject
        (
            "sGamma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("one", dimless, 0)
    );

    const surfaceScalarField& magSf = mesh.magSf();
    const scalarField magFaceAreas = mag(mesh.faceAreas());

    sGamma.internalField() =
        magSf.internalField()/
        scalarField::subField(magFaceAreas, mesh.nInternalFaces());

    forAll (mesh.boundary(), patchI)
    {
        if (!isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
        {
            sGamma.boundaryField()[patchI] =
                magSf.boundaryField()[patchI]/
                mesh.boundary()[patchI].patchSlice(magFaceAreas);

            gamma.boundaryField()[patchI] =
                sGamma.boundaryField()[patchI];
        }
    }

    sGamma.write();
    gamma.write();


    // Create dead cells set
    {
        cellSet
        (
            mesh,
            "deadCells",
            deadCellsHash
        ).write();
    }
    
    // Check consistency of face area vectors

    Info<< nl << "Calculating divSf" << endl;
    volVectorField divSf
    (
        "divSf",
        fvc::surfaceIntegrate(mesh.Sf())
    );
    divSf.write();

    // Check divergence of face area vectors. Note: scale by the volume
    // to avoid bias towards small cells.  HJ, 13/Mar/2019
    scalarField magDivSf = mag(divSf)().internalField()*mesh.V().field();

    Info<< "Face areas divergence (min, max, average): "
        << "(" << min(magDivSf) << " " << max(magDivSf)
        << " " << average(magDivSf) << ")"
        << endl;

    if (max(magDivSf) > primitiveMesh::closedThreshold_)
    {
        WarningIn("writeIbMasks")
            << "Possible problem with immersed boundary face area vectors: "
            << max(magDivSf)
            << endl;

        scalar maxOpenCell = 0;
        label maxOpenCellIndex = -1;

        forAll (magDivSf, cellI)
        {
            if (magDivSf[cellI] > maxOpenCell)
            {
                maxOpenCell = magDivSf[cellI];
                maxOpenCellIndex = cellI;
            }

            if (magDivSf[cellI] > 1e-9)
            {
                Info<< "Open cell " << cellI << ": " << magDivSf[cellI]
                    << " gamma: " << gamma[cellI] << endl;
            }
        }

        const surfaceVectorField& Sf = mesh.Sf();

        const labelList& openCellFaces = mesh.cells()[maxOpenCellIndex];

        scalarField openCellFaceGamma(openCellFaces.size(), scalar(-1));

        vectorField openFaceAreas
        (
            IndirectList<vector>(mesh.faceAreas(), openCellFaces)()
        );

        vectorField adjustedFaceAreas(openCellFaces.size());

        forAll (openCellFaces, cfI)
        {
            const label& faceI = openCellFaces[cfI];

            if (mesh.isInternalFace(faceI))
            {
                openCellFaceGamma[cfI] = sGamma.internalField()[faceI];

                adjustedFaceAreas[cfI] = Sf.internalField()[faceI];
            }
            else
            {
                const label patchI = mesh.boundaryMesh().whichPatch(faceI);

                const label patchFaceI =
                    mesh.boundaryMesh()[patchI].whichFace(faceI);

                openCellFaceGamma[cfI] =
                    sGamma.boundaryField()[patchI][patchFaceI];

                adjustedFaceAreas[cfI] = Sf.boundaryField()[patchI][patchFaceI];
            }
        }

        // Find faces on IB patches
        vectorField ibVectors(mesh.boundary().size());
        label nIbVectors = 0;

        forAll (mesh.boundary(), patchI)
        {
            if (isA<immersedBoundaryFvPatch>(mesh.boundary()[patchI]))
            {
                const labelList& ibpFC = mesh.boundary()[patchI].faceCells();

                forAll (ibpFC, ibpFCI)
                {
                    if (ibpFC[ibpFCI] == maxOpenCellIndex)
                    {
                        ibVectors[nIbVectors] =
                            mesh.Sf().boundaryField()[patchI][ibpFCI];

                        nIbVectors++;
                    }
                }
            }
        }

        ibVectors.setSize(nIbVectors);

        Pout<< "Max open cell index: " << maxOpenCellIndex
            << " magDivSf = " << maxOpenCell << nl
            << "faces: " << openCellFaces << nl
            << " original areas: " << openFaceAreas << nl
            << "sGamma: " << openCellFaceGamma << nl
            << "adjusted areas: " << adjustedFaceAreas << nl
            << "cut face areas: " << ibVectors << nl
            << "Sum normal areas: " << sum(openFaceAreas) << nl
            << "Sum iB areas: " << sum(ibVectors) << nl
            << endl;
    }

    Info<< endl;
}


// ************************************************************************* //
