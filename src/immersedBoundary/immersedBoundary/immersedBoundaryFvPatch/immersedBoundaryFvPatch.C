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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatch.H"
#include "foamTime.H"
#include "fvBoundaryMesh.H"
#include "fvMesh.H"
#include "SortableList.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "treeBoundBox.H"
#include "treeDataCell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(immersedBoundaryFvPatch, 0);

    addToRunTimeSelectionTable(fvPatch, immersedBoundaryFvPatch, polyPatch);
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::radiusFactor_
(
    "immersedBoundaryRadiusFactor",
    3.5
);


const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::angleFactor_
(
    "immersedBoundaryAngleFactor",
    80
);


const Foam::debug::optimisationSwitch
Foam::immersedBoundaryFvPatch::maxCellCellRows_
(
    "immersedBoundaryMaxCellCellRows",
    4
);


const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::distFactor_
(
    "immersedBoundaryDistFactor",
    1.5
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::immersedBoundaryFvPatch::csEst() const
{
    return Foam::max(100, 0.1*mesh_.nCells());
}


void Foam::immersedBoundaryFvPatch::makeGamma() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeGamma() const")
            << "creating fluid cells indicator "
            << "for immersed boundary" << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeGamma() const")
            << "fluid cells indicator already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    gammaPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGamma" + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1)
        );

    scalarField& gammaI = gammaPtr_->internalField();

    // Start from gammaExt and mark IB cells and inactive
    gammaI = gammaExt().internalField();

    const labelList& ibc = ibCells();

    // Remove IB cells from gammaExt to create gamma
    forAll (ibc, cellI)
    {
        gammaI[ibc[cellI]] = 0;
    }

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones
    gammaPtr_->boundaryField().evaluateCoupled();

    forAll (gammaPtr_->boundaryField(), patchI)
    {
        if (!gammaPtr_->boundaryField()[patchI].coupled())
        {
            gammaPtr_->boundaryField()[patchI] =
                gammaPtr_->boundaryField()[patchI].patchInternalField();
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeGammaExt() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeGammaExt() const")
            << "creating extended fluid cells indicator "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (gammaExtPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeGammaExt() const")
            << "extended fluid cells indicator already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    gammaExtPtr_ =
        new volScalarField
        (
            IOobject
            (
                "ibGammaExt"  + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("1", dimless, 1)
        );

    scalarField& gammaExtI = gammaExtPtr_->internalField();

    const vectorField& C = mesh_.cellCentres();

    // Mark cells that are inside or outside of the triangular surface
    boolList inside = ibPolyPatch_.triSurfSearch().calcInside(C);

    // Adjust selection of cells: inside or outside of immersed boundary
    if (internalFlow())
    {
        Info<< "Internal flow" << endl;
        forAll (gammaExtI, cellI)
        {
            if (!inside[cellI])
            {
                gammaExtI[cellI] = 0;
            }
        }
    }
    else
    {
        Info<< "External flow" << endl;
        forAll (gammaExtI, cellI)
        {
            if (inside[cellI])
            {
                gammaExtI[cellI] = 0;
            }
        }
    }

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones
    gammaExtPtr_->boundaryField().evaluateCoupled();
}


void Foam::immersedBoundaryFvPatch::makeSGamma() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeSGamma() const")
            << "creating fluid faces indicator "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (sGammaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeSGamma() const")
            << "fluid faces indicator already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Note: change in algorithm
    // 1) First, mark all faces as dead
    // 2) Mark faces as live if they are between two live cells
    sGammaPtr_ =
        new surfaceScalarField
        (
            IOobject
            (
                "sGamma"  + name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimless, 0),
            calculatedFvsPatchScalarField::typeName
        );

    // Get access to components of sGamma
    scalarField& sGammaI = sGammaPtr_->internalField();

    surfaceScalarField::GeometricBoundaryField& sGammaPatches =
        sGammaPtr_->boundaryField();

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    // Live cells indicator
    const volScalarField& g = gamma();

    // Extended live cells indicator
    const volScalarField& gExt = gammaExt();
    const scalarField& gExtIn = gExt.internalField();

    const volScalarField::GeometricBoundaryField& gExtPatches =
        gExt.boundaryField();

    volScalarField gIb = gExt - g;
    const scalarField& gIbIn = gIb.internalField();

    const volScalarField::GeometricBoundaryField& gIbPatches =
        gIb.boundaryField();

    // Internal faces: flux is live between all active and IB cells
    forAll (sGammaI, faceI)
    {
        // If both cells are live, flux is live
        // Note: checking 0:1 indicator
        if
        (
            gExtIn[owner[faceI]] > SMALL
         && gExtIn[neighbour[faceI]] > SMALL
        )
        {
            sGammaI[faceI] = 1;
        }
    }

    // Kill fluxes between two IB cells
    forAll (sGammaI, faceI)
    {
        // If both cells are live, flux is live
        // Note: checking 0:1 indicator
        if
        (
            gIbIn[owner[faceI]] > SMALL
         && gIbIn[neighbour[faceI]] > SMALL
        )
        {
            sGammaI[faceI] = 0;
        }
    }

    forAll (gExtPatches, patchI)
    {
        if (gExtPatches[patchI].coupled())
        {
            scalarField& gP = sGammaPatches[patchI];

            // For coupled patches, check gammaExt
            scalarField gammaOwn = gExtPatches[patchI].patchInternalField();

            scalarField gammaNei = gExtPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                // Note: checking 0:1 indicator
                if
                (
                    gammaOwn[faceI] > SMALL
                 && gammaNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 1;
                }
            }

            // For coupled patches, kill IB
            scalarField gammaIbOwn = gIbPatches[patchI].patchInternalField();

            scalarField gammaIbNei = gIbPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                // Note: checking 0:1 indicator
                if
                (
                    gammaIbOwn[faceI] > SMALL
                 && gammaIbNei[faceI] > SMALL
                )
                {
                    gP[faceI] = 0;
                }
            }
        }
        else
        {
            // For regular patches, check live cells only to achieve
            // correct global mass adjustment.
            // HJ, 21/May/2012
            scalarField gammaFc =
                g.boundaryField()[patchI].patchInternalField();

            scalarField& gP = sGammaPatches[patchI];

            forAll (gammaFc, faceI)
            {
                // Note: checking 0:1 indicator
                if (gammaFc[faceI] > SMALL)
                {
                   gP[faceI] = 1;
                }
            }
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeIbCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCells() const")
            << "create list of cells next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCells() const")
            << "list of cells next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Collect IB cells in hashSet to eliminate duplicates
    labelHashSet ibCellSet(csEst());

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    // IB cell is identified by having a dead neighbour
    forAll (neighbour, faceI)
    {
        // Note: checking 0:1 indicator
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            // Note: checking 0:1 indicator
            if (gammaExtI[owner[faceI]] > SMALL)
            {
                // Owner is live: IB cell
                // Search not needed: duplicates automatically filtered
                // HJ, 2/May/2017
                ibCellSet.insert(owner[faceI]);
            }
            else
            {
                // Neighbour is live: IB cell
                // Search not needed: duplicates automatically filtered
                // HJ, 2/May/2017
                ibCellSet.insert(neighbour[faceI]);
            }
        }
    }

    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gE.boundaryField()[patchI].patchInternalField();

            // Get gamma from other side of a coupled patch
            scalarField gammaExtNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells =
                mesh_.boundary()[patchI].faceCells();

            forAll (gammaExtOwn, faceI)
            {
                // Note: checking 0:1 indicator
                if
                (
                    mag(gammaExtNei[faceI] - gammaExtOwn[faceI])
                  > SMALL
                )
                {
                    if (gammaExtOwn[faceI] > SMALL)
                    {
                        // Search not needed: duplicates automatically filtered
                        // HJ, 2/May/2017
                        ibCellSet.insert(fCells[faceI]);
                    }
                    // Bugfix: Removed special handling for cyclic patch
                }
            }
        }
    }

    ibCellsPtr_ = new labelList(ibCellSet.sortedToc());

    Pout<< "Number of IB cells: " << ibCellsPtr_->size() << endl;
}


void Foam::immersedBoundaryFvPatch::addIbCornerCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::addIbCornerCells() const")
            << "add cells next to sharp corner "
            << "for immersed boundary " << name()
            << endl;
    }

    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const cellList& meshCells = mesh_.cells();
    const unallocLabelList& own = mesh_.owner();
    const unallocLabelList& nei = mesh_.neighbour();

    label nCornerCells = 0;

    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();

    labelList cornerCells;

    do
    {
        // Access to derived data from do-loop: this will be re-calculated
        // HJ, 21/May/2012
        const labelList& ibc = ibCells();

        // Note: the algorithm is originally written with inward-facing normals
        // and subsequently changed: IB surface normals point outwards
        // HJ, 21/May/2012
        const vectorField& ibn = ibNormals();

        const scalarField& gammaI = gamma().internalField();

        labelHashSet cornerIbCellSet(csEst());

        const labelList& ibf = ibFaces();

        forAll (ibf, faceI)
        {
            const label& ownCell = own[ibf[faceI]];
            const label& neiCell = nei[ibf[faceI]];

            label liveCell = -1;

            // Note: checking 0:1 indicator
            if (gammaI[ownCell] > SMALL)
            {
                liveCell = ownCell;
            }
            else
            {
                liveCell = neiCell;
            }

            scalar delta = cellSize(liveCell);
            vector span(2*delta, 2*delta, 2*delta);

            pointIndexHit pih = tss.nearest(C[liveCell], span);

            if (pih.hit())
            {
                vector n =
                    triSurfaceTools::surfaceNormal
                    (
                        ibPolyPatch_.ibMesh(),
                        pih.index(),
                        pih.hitPoint()
                    );

                scalar totalArea = 0;
                {
                    const labelList& cellFaces = meshCells[liveCell];

                    forAll (cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];

                        vector curSf = Sf[curFace];

                        if ((curSf & (Cf[curFace] - C[liveCell])) < 0)
                        {
                            curSf *= -1;
                        }

                        if ((curSf & n) > 0)
                        {
                            totalArea += (curSf & n);
                        }
                    }
                }

                scalar area = 0;
                {
                    const labelList& cellFaces = meshCells[liveCell];

                    forAll (cellFaces, faceI)
                    {
                        label curFace = cellFaces[faceI];

                        vector curSf = Sf[curFace];

                        label neiCell = -1;

                        if (mesh_.isInternalFace(curFace))
                        {
                            if (own[curFace] == liveCell)
                            {
                                neiCell = nei[curFace];
                            }
                            else
                            {
                                neiCell = own[curFace];
                            }
                        }

                        label ibCell = findIndex(ibc, neiCell);

                        if (ibCell != -1)
                        {
                            // Note that ibn points outwards
                            if ((-ibn[ibCell] & n) > 0)
                            {
                                area += mag(curSf & n);
                            }
                        }
                    }
                }

                if (area/totalArea < 0.5)
                {
                    // Insert corner cell without checking
                    cornerIbCellSet.insert(liveCell);
                }
            }
        }

        cornerCells = cornerIbCellSet.toc();
        ibCellsPtr_->append(cornerCells);
        nCornerCells += cornerCells.size();

        deleteDemandDrivenData(gammaPtr_);

        deleteDemandDrivenData(ibFacesPtr_);
        deleteDemandDrivenData(ibFaceCellsPtr_);
        deleteDemandDrivenData(ibFaceFlipsPtr_);

        deleteDemandDrivenData(ibPointsPtr_);
        deleteDemandDrivenData(ibNormalsPtr_);
        deleteDemandDrivenData(hitFacesPtr_);
        deleteDemandDrivenData(ibSamplingPointsPtr_);

        deleteDemandDrivenData(ibSamplingWeightsPtr_);
        deleteDemandDrivenData(ibSamplingProcWeightsPtr_);
    }
    while (cornerCells.size() > 0);
}


void Foam::immersedBoundaryFvPatch::makeIbFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbFaces() const")
            << "create list of faces next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibFacesPtr_ || ibFaceCellsPtr_ || ibFaceFlipsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbFaces() const")
            << "list of faces next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Mark IB cells with their index
    const labelList& ibc = ibCells();

    labelList ibCellIndicator(mesh_.nCells(), -1);

    forAll (ibc, ibcI)
    {
        ibCellIndicator[ibc[ibcI]] = ibcI;
    }

    dynamicLabelList ibF(2*ibc.size());
    dynamicLabelList ibFC(2*ibc.size());
    DynamicList<bool> ibFF(2*ibc.size());

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            ibCellIndicator[owner[faceI]] == -1
         && ibCellIndicator[neighbour[faceI]] > -1
        )
        {
            // Owner is live, neighbour IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[neighbour[faceI]]);
            ibFF.append(false);
        }
        else if
        (
            ibCellIndicator[owner[faceI]] > -1
         && ibCellIndicator[neighbour[faceI]] == -1
        )
        {
            // Neighbour is live, owner IB.  Its IB index is in
            // ibCellIndicator
            ibF.append(faceI);
            ibFC.append(ibCellIndicator[owner[faceI]]);
            ibFF.append(true);
        }
    }

    const volScalarField::GeometricBoundaryField& gammaPatches =
        gamma().boundaryField();

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh_.boundary()[patchI].faceCells();
        const label start = mesh_.boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {
            scalarField gammaOwn =
                gammaPatches[patchI].patchInternalField();

            scalarField gammaNei =
                gammaPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, patchFaceI)
            {
                if
                (
                    mag(gammaNei[patchFaceI] - gammaOwn[patchFaceI]) > SMALL
                )
                {
                    if (ibCellIndicator[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(ibCellIndicator[fc[patchFaceI]]);
                        ibFF.append(false);
                    }
                    else
                    {
                        // Neighbour cell is IB
                        ibF.append(start + patchFaceI);
                        ibFC.append(-1);
                        ibFF.append(true);
                    }
                }
            }
        }
    }

    // Pack the data
    ibF.shrink();
    ibFC.shrink();
    ibFF.shrink();

    ibFacesPtr_ = new labelList(ibF.xfer());

    ibFaceCellsPtr_ = new labelList(ibFC.xfer());
    ibFaceFlipsPtr_ = new boolList(ibFF.xfer());
}


void Foam::immersedBoundaryFvPatch::makeIbInsideFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbInsideFaces() const")
            << "create list of faces next to immersed boundary "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibInsideFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbInsideFaces() const")
            << "list of faces next to immersed boundary already exist "
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Create face hash set with estimated size
    labelHashSet ibInsideFSet(Foam::max(100, 0.1*mesh_.nFaces()));

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    forAll (neighbour, faceI)
    {
        // Note: checking 0:1 indicator
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            ibInsideFSet.insert(faceI);
        }
    }

    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gE.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            label patchStart = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                // Note: checking 0:1 indicator
                if
                (
                    mag(gammaNei[faceI] - gammaOwn[faceI]) > SMALL
                )
                {
                    // Search not needed: duplicates automatically filtered
                    // HJ, 2/May/2017
                    ibInsideFSet.insert(patchStart + faceI);
                }
            }
        }
    }

    ibInsideFacesPtr_ = new labelList(ibInsideFSet.sortedToc());
}


void Foam::immersedBoundaryFvPatch::makeIbInternalFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbInternalFaces() const")
            << "create list of faces next to immersed boundary"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibInternalFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbInternalFaces() const")
            << "list of faces next to immersed boundary already exist"
            << abort(FatalError);
    }

    // Create face hash set with estimated size
    labelHashSet ibInternalFacesSet(Foam::max(100, 0.1*mesh_.nFaces()));

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    volScalarField gammaTmp
    (
        "gammaTmp",
        gammaExt() - gamma()
    );
    const scalarField& gammaTmpI = gammaTmp.internalField();

    forAll (neighbour, faceI)
    {
        // Note: checking 0:1 indicator
        if
        (
            (gammaTmpI[neighbour[faceI]] > SMALL)
         && (gammaTmpI[owner[faceI]] > SMALL)
        )
        {
            ibInternalFacesSet.insert(faceI);
        }
    }

    forAll (gammaTmp.boundaryField(), patchI)
    {
        if (gammaTmp.boundaryField()[patchI].coupled())
        {
            scalarField gammaOwn =
                gammaTmp.boundaryField()[patchI].patchInternalField();

            scalarField gammaNei =
                gammaTmp.boundaryField()[patchI].patchNeighbourField();

            label patchStart = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                // Note: checking 0:1 indicator
                if
                (
                    (gammaNei[faceI] > SMALL)
                 && (gammaOwn[faceI] > SMALL)
                )
                {
                    // Search not needed: duplicates automatically filtered
                    // HJ, 2/May/2017
                    ibInternalFacesSet.insert(patchStart + faceI);
                }
            }
        }
    }

    ibInternalFacesPtr_ = new labelList(ibInternalFacesSet.sortedToc());
}


void Foam::immersedBoundaryFvPatch::makeIbPointsAndNormals() const
{
    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
        )   << "create immersed  boundary points and normals "
            << "for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibPointsPtr_ || ibNormalsPtr_ || hitFacesPtr_ || ibSamplingPointsPtr_)
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
        )
            << "immersed boundary points and normals already exist"
            << "for immersed boundary " << name()
            << abort(FatalError);
    }

    // Find representative cell dimension for IB cells
    const labelList& ibc = ibCells();

    scalarField delta(ibc.size());

    forAll (delta, cellI)
    {
        delta[cellI] = cellSize(ibc[cellI]);
    }

    // Find nearest triSurface point for each interface cell centre

    // Allocate storage
    ibPointsPtr_ = new vectorField(ibc.size(), vector::zero);
    vectorField& ibPoints = *ibPointsPtr_;

    ibNormalsPtr_ = new vectorField(ibc.size(), vector::zero);
    vectorField& ibNormals = *ibNormalsPtr_;

    hitFacesPtr_ = new labelList(ibc.size(), -1);
    labelList& ibHitFaces = *hitFacesPtr_;

    ibSamplingPointsPtr_ = new vectorField(ibc.size(), vector::zero);
    vectorField& ibSamplingPoints = *ibSamplingPointsPtr_;

    // Get IB cell centres by subsetting
    vectorField ibCellCentres(mesh_.cellCentres(), ibc);

    // Get surface search
    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();

    forAll (ibc, cellI)
    {
        // Adjust search span if needed.  HJ, 14/Dec/2012
        vector span
        (
            2*radiusFactor_()*delta[cellI],
            2*radiusFactor_()*delta[cellI],
            2*radiusFactor_()*delta[cellI]
        );

        pointIndexHit pih = tss.nearest(ibCellCentres[cellI], span);

        if (pih.hit())
        {
            ibPoints[cellI] = pih.hitPoint();
            ibNormals[cellI] =
                triSurfaceTools::surfaceNormal
                (
                    ibPolyPatch_.ibMesh(),
                    pih.index(),
                    pih.hitPoint()
                );

            // Note: ibNormals point OUT of the domain
            if (!internalFlow())
            {
                ibNormals[cellI] *= -1;
            }

            ibHitFaces[cellI] = pih.index();
        }
        else
        {
            FatalErrorIn
            (
                "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
            )   << "Can't find nearest triSurface point for cell "
                << ibc[cellI] << ", "
                << mesh_.cellCentres()[ibc[cellI]]
                << ".  Hit data = " << pih << nl
                << abort(FatalError);
        }

        if
        (
            mesh_.nGeometricD() < 3
         && mag(ibCellCentres[cellI].z() - ibPoints[cellI].z()) > SMALL
        )
        {
            WarningIn
            (
                "immersedBoundaryFvPatch::makeIbPointsAndNormals() const"
            )   << "Intersection point is not on symmetry plane " << nl
                << "C = " << ibCellCentres[cellI]
                <<  " D = " <<  ibPoints[cellI] << nl
                << "for 2-D geometry.  Adjusting" << endl;

               ibPoints[cellI].z() = ibCellCentres[cellI].z();
        }
    }

    // Calculate sampling points locations
    ibSamplingPoints = ibPoints + distFactor_()*(ibCellCentres - ibPoints);
}


void Foam::immersedBoundaryFvPatch::makeIbCellCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCellCells() const")
            << "create neighbour cells for ib cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if
    (
        ibCellCellsPtr_
     || ibProcCellsPtr_
     || ibMapPtr_
     || ibProcCentresPtr_
     || ibProcGammaPtr_
     || ibCellProcCellsPtr_
    )
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbCellCells() const"
        )   << "cell-cell addressing already exists"
            << abort(FatalError);
    }

    // Get mesh cells
    const cellList& meshCells = mesh_.cells();

    // In geometry initialisation, fields are not available: use raw mesh data
    // HJ after ZT, 6/Dec/2012
    const vectorField& C = mesh_.cellCentres();

    // Get IB cells
    const labelList& ibc = ibCells();

    // Allocate storage

    // Local IB support: cellCells
    ibCellCellsPtr_ = new labelListList(ibc.size());
    labelListList& cellCells = *ibCellCellsPtr_;

    // Dummy map: for parallel run it will be reset below
    {
        labelListList zeroList;

        ibMapPtr_ = new mapDistribute
        (
            0,
            zeroList,
            zeroList,
            false            // Do not reuse
        );
    }

    // Cells that support interpolation on other
    ibProcCellsPtr_ = new labelListList(Pstream::nProcs());
    labelListList& procCells = *ibProcCellsPtr_;

    // Centres of other supporting processors.  Sized only in parallel run
    ibProcCentresPtr_ = new vectorList();
    vectorList& procCentres = *ibProcCentresPtr_;

    // Gamma from other processorss.  Sized only in parallel run
    ibProcGammaPtr_ = new scalarList();
    scalarList& procGamma = *ibProcGammaPtr_;

    // Neighbour IB cell-proc-cell addressing
    ibCellProcCellsPtr_ = new labelListList(ibc.size());
    labelListList& cellProcCells = *ibCellProcCellsPtr_;


    // Get characteristic cell size for IB cells
    const scalarField& cellSizes = ibCellSizes();

    // Calculate angle limit
    scalar angleLimit = -Foam::cos(angleFactor_()*mathematicalConstant::pi/180);

    // Calculate radius for all ibCells
    scalarField rM = radiusFactor_()*cellSizes;

    // Get IB points
    const vectorField& ibp = ibPoints();

    // Note: the algorithm is originally written with inward-facing normals
    // and subsequently changed: IB surface normals point outwards
    // HJ, 21/May/2012
    const vectorField& ibn = ibNormals();

    forAll (ibc, cellI)
    {
        // Get IB cell centre
        const vector curIbCentre = C[ibc[cellI]];

        // Collect extended neighbourhood for search
        labelList curCells;

        findCellCells
        (
            curIbCentre,
            ibc[cellI],
            curCells
        );

        // Filter the cellCells for the direction angle
        cellCells[cellI] = labelList(curCells.size(), -1);

        label cI = 0;

        forAll (curCells, i)
        {
            label curCell = curCells[i];

            // Collect the cells within rM of the fitting cell
            if (mag(C[curCell] - C[ibc[cellI]]) <= rM[cellI])
            {
                // Within distance.  Search direction
                vector dir = (C[curCell] - ibp[cellI]);
                dir /= mag(dir) + SMALL;

                // Change of sign of normal.  HJ, 21/May/2012
                if ((-ibn[cellI] & dir) >= angleLimit)
                {
                    cellCells[cellI][cI++] = curCell;
                }
            }
        }

        cellCells[cellI].setSize(cI);
    }

    // Establish immersed boundary support spanning across processor
    // boundaries
    // Algorithm:
    // 1) Find all local IB clusters that touch a processor boundary.
    //    For each such cluster collect the cell index, centroid and rMax
    // 2) Using gather-scatter distribute the search data to all processors
    // 3) For each processor, perform a local search for all other processors
    //    (but not for self)

    // Check the need for parallel communication
    // If the cellCell cluster touches a processor boundary, it may be possible
    // to locate interpolation support on other processors
    if (Pstream::parRun())
    {
        // Set the fields that exist only for a parallel run.  They were
        // created empty
        procCells.setSize(Pstream::nProcs());

        // Find immersed boundary cells next to processor boundaries
        labelHashSet procIbCellsSet(csEst());

        forAll (ibc, cellI)
        {
            const labelList& curCellCells = cellCells[cellI];

            // If there are local neighbours, check them
            if (!curCellCells.empty())
            {
                forAll (curCellCells, cI)
                {
                    const labelList& faces = meshCells[curCellCells[cI]];

                    bool foundProcessorFace = false;

                    forAll (faces, faceI)
                    {
                        if (!mesh_.isInternalFace(faces[faceI]))
                        {
                            // Face is not internal.  Find patch
                            label patchID =
                                mesh_.boundaryMesh().whichPatch(faces[faceI]);

                            // If a processor patch is hit, prepare search
                            // on other processors
                            // HJ, Rewrite: this requires generalisation
                            // HJ, 1/Oct/2017
                            if
                            (
                                isA<processorPolyPatch>
                                (
                                    mesh_.boundaryMesh()[patchID]
                                )
                            )
                            {
                                foundProcessorFace = true;
                            }
                        }
                    }

                    if (foundProcessorFace)
                    {
                        procIbCellsSet.insert(cellI);
                        break;
                    }
                }
            }
            else
            {
                // No local neighbours detected.  Check IB cell
                const labelList& faces = meshCells[ibc[cellI]];

                bool foundProcessorFace = false;

                forAll (faces, faceI)
                {
                    if (!mesh_.isInternalFace(faces[faceI]))
                    {
                        // Face is not internal.  Find patch
                        label patchID =
                            mesh_.boundaryMesh().whichPatch(faces[faceI]);

                        // If a processor patch is hit, prepare search
                        // on other processors
                        // HJ, Rewrite: this requires generalisation
                        // HJ, 1/Oct/2017
                        if
                        (
                            isA<processorPolyPatch>
                            (
                                mesh_.boundaryMesh()[patchID]
                            )
                        )
                        {
                            foundProcessorFace = true;
                        }
                    }
                }

                if (foundProcessorFace)
                {
                    procIbCellsSet.insert(cellI);
                }
            }
        }

        // Collect all IB cells that require processor search
        labelList procIbCells = procIbCellsSet.sortedToc();

        // Note: new gather-scatter operations
        // HJ, 11/Aug/2016

        // Note: possible optimisation
        // It is possible to avoid sending all IB cell seeds to all processors
        // if the seed point is "sufficiently far" from the processor mesh
        // bounding box.
        // This reduces communication, but does not affect search time
        // For details, see equivalent overset mesh code in oversetRegion.C
        // HJ, 4/Oct/2017

        // Send and receive ibc centres and radii
        // This is used for search on other processors
        vectorListList ctrs(Pstream::nProcs());

        ctrs[Pstream::myProcNo()].setSize(procIbCells.size());
        vectorList& centres = ctrs[Pstream::myProcNo()];

        forAll (centres, cellI)
        {
            centres[cellI] = C[ibc[procIbCells[cellI]]];
        }

        // Broadcast ibc centres
        Pstream::gatherList(ctrs);
        Pstream::scatterList(ctrs);

        // Send and receive IB procRMax
        scalarListList procRMax(Pstream::nProcs());
        procRMax[Pstream::myProcNo()] = scalarField(rM, procIbCells);

        // Broadcast procRMax
        // This is used for search on other processors
        Pstream::gatherList(procRMax);
        Pstream::scatterList(procRMax);

        // Send and receive IB procIbn
        vectorListList procIbn(Pstream::nProcs());
        procIbn[Pstream::myProcNo()] = vectorField(ibn, procIbCells);

        // Broadcast procRMax
        // This is used for search on other processors
        Pstream::gatherList(procIbn);
        Pstream::scatterList(procIbn);

        // Get search tree
        const indexedOctree<treeDataCell>& tree = cellSearch();

        const scalar span = tree.bb().mag();

        // Using octree search, establish send-receive maps for parallel data
        // collection.  HJ, 4/Oct/2017

        // It is possible that an octree is empty.  If there are no eligible
        // cells on this processor, do not search
        if (!tree.nodes().empty())
        {
            // Search all processor clusters from other processors
            // to find local IB support
            for (label procI = 0; procI < Pstream::nProcs(); procI++)
            {
                // Search all processor apart from self
                if (procI != Pstream::myProcNo())
                {
                    labelHashSet procCellSet(csEst());

                    // Get points to search
                    const vectorList& curCtrs = ctrs[procI];
                    const scalarList& curRMax = procRMax[procI];

                    forAll (curCtrs, cellI)
                    {
                        // Get current point
                        const point& curP = curCtrs[cellI];

                        // Find nearest cell with octree. Note: octree only
                        // contains eligible cells.  HJ, 10/Jan/2015.
                        const pointIndexHit pih = tree.findNearest(curP, span);

                        if (pih.hit())
                        {
                            // Get index obtained by octree
                            const label nearestCellID = pih.index();

                            // Valid hit found.  Check radius

                            // Calculate radius
                            scalar R = mag(C[nearestCellID] - curP);

                            if (R < curRMax[cellI])
                            {
                                // Insert cell as support
                                // Search not needed: duplicates
                                // automatically filtered
                                // HJ, 2/May/2017
                                procCellSet.insert(nearestCellID);

                                // Search neighbourhood
                                labelList tmpCellList;

                                // Collect extended neighbourhood
                                // for search
                                findCellCells
                                (
                                    curP,
                                    nearestCellID,
                                    tmpCellList
                                );

                                forAll (tmpCellList, cI)
                                {
                                    const scalar r =
                                        mag
                                        (
                                            C[tmpCellList[cI]]
                                          - curP
                                        );

                                    if (r <= procRMax[procI][cellI])
                                    {
                                        // Within distance.
                                        // Search direction
                                        vector dir = (C[nearestCellID] - curP);
                                        dir /= mag(dir) + SMALL;

                                        // Change of sign of normal.
                                        // HJ, 21/May/2012
                                        if
                                        (
                                            (-procIbn[procI][cellI] & dir)
                                         >= angleLimit
                                        )
                                        {
                                            // Search not needed: duplicates
                                            // automatically filtered
                                            // HJ, 2/May/2017
                                            procCellSet.insert
                                            (
                                                tmpCellList[cI]
                                            );
                                        }
                                    }
                                }
                            }
                        }
                    }

                    // Collected all local cells that may act as support for
                    // other processor IB points
                    procCells[procI] = procCellSet.toc();
                }
            }
        }

        // Note:
        // procCells is the send map.  The receive map is calculated by
        // finding out how many cells are passed from each processor
        // and ordering them in processor order
        // This will be a square matrix in which each processor fills its row
        labelListList nAllProcUsed(Pstream::nProcs());
        labelList& curProcUsed = nAllProcUsed[Pstream::myProcNo()];
        curProcUsed.setSize(Pstream::nProcs(), 0);

        forAll (curProcUsed, procI)
        {
            curProcUsed[procI] = procCells[procI].size();
        }

        // Gather/scatter number of used points going to each processor from
        // each processor so that all processors have all necessary information
        // when creating the map distribute tool for distributing points
        Pstream::gatherList(nAllProcUsed);
        Pstream::scatterList(nAllProcUsed);

        // Create components for mapDistribute
        // Numbering of received points: all points from proc 1,
        // followed by all points from proc 2 etc.

        // Record how to assemble points from which processor
        labelListList constructMap(Pstream::nProcs());

        label constructSize = 0;

        forAll (nAllProcUsed, procI)
        {
            // Get column size
            const label colSize = nAllProcUsed[procI][Pstream::myProcNo()];

            // Fill the map
            labelList& curConstructMap = constructMap[procI];

            curConstructMap.setSize(colSize);

            forAll (curConstructMap, i)
            {
                curConstructMap[i] = constructSize;
                constructSize++;
            }
        }

        // Create mapDistribute
        deleteDemandDrivenData(ibMapPtr_);

        ibMapPtr_ = new mapDistribute
        (
            constructSize,
            procCells,
            constructMap,
            false            // Do not reuse
        );
        const mapDistribute& ibm = *ibMapPtr_;

        // Use mapDistribute to assemble other data

        // Centres from other processors
        procCentres = C;
        ibm.distribute(procCentres);

        // Gamma from other processors
        procGamma = gamma().internalField();
        ibm.distribute(procGamma);

        // Algorithm
        // At this stage, all possible cells providing support for IB points
        // on other processor have been collected and broadcast.
        // Go through all ib cells and select the points
        // that are used for individual local IB cells.

        forAll (ibc, cellI)
        {
            // Get IB cell centre
            const vector curIbCentre = C[ibc[cellI]];

            labelList& curCellProcCells = cellProcCells[cellI];

            curCellProcCells.setSize(procCentres.size());

            // Count number of proceCells to use
            label cI = 0;

            // Check all procCentres
            // Note
            // This is a problematic algorithm because it is squared in
            // the number of IB cells and procCells, neither of which
            // can be controlled
            // Formally, an IB cell which is further away from the processor
            // boundary may also have support on other processors.
            // To handle this, all IB cells are checked against all other
            // procCells that support

            forAll (procCentres, i)
            {
                // Collect the cells within rM of the fitting cell
                if (mag(procCentres[i] - curIbCentre) <= rM[cellI])
                {
                    // Within distance.  Search direction
                    vector dir = (procCentres[i] -  ibp[cellI]);
                    dir /= mag(dir) + SMALL;

                    // Change of sign of normal.  HJ, 21/May/2012
                    if ((-ibn[cellI] & dir) >= angleLimit)
                    {
                        curCellProcCells[cI++] = i;
                    }
                }
            }

            // Resize the cell list
            curCellProcCells.setSize(cI);
        }
    }

    // if (debug)
    {
        // Check size of cellCells, cellProcCells and sum for all ibCells
        label minCC = INT_MAX;
        label minProcCC = INT_MAX;
        label minTotCC = INT_MAX;

        label maxCC = 0;
        label maxProcCC = 0;
        label maxTotCC = 0;

        forAll (cellCells, ibpI)
        {
            minCC = Foam::min(minCC, cellCells[ibpI].size());
            minProcCC = Foam::min(minProcCC, cellProcCells[ibpI].size());
            minTotCC = Foam::min
            (
                minTotCC,
                cellCells[ibpI].size() + cellProcCells[ibpI].size()
            );

            maxCC = Foam::max(maxCC, cellCells[ibpI].size());
            maxProcCC = Foam::max(maxProcCC, cellProcCells[ibpI].size());
            maxTotCC = Foam::max
            (
                maxTotCC,
                cellCells[ibpI].size() + cellProcCells[ibpI].size()
            );
        }

        reduce(minCC, minOp<label>());
        reduce(minProcCC, minOp<label>());
        reduce(minTotCC, minOp<label>());

        reduce(maxCC, maxOp<label>());
        reduce(maxProcCC, maxOp<label>());
        reduce(maxTotCC, maxOp<label>());

        Info<< "Min/Max support CC, procCC, total = ("
            << minCC << token::SPACE << minProcCC << token::SPACE << minTotCC
            << ", "
            << maxCC << token::SPACE << maxProcCC << token::SPACE << maxTotCC
            << ")" << endl;
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadCells() const")
            << "create list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadCells() const")
            << "list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaExtI = gammaExt().internalField();

    deadCellsPtr_ = new labelList(label(sum(scalar(1) - gammaExtI)), -1);
    labelList& deadCells = *deadCellsPtr_;

    label counter = 0;

    forAll (gammaExtI, cellI)
    {
        if (gammaExtI[cellI] < SMALL)
        {
            deadCells[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadCellsExt() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadCellsExt() const")
            << "create extended list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadCellsExtPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadCellsExt() const")
            << "extended list of dead cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    deadCellsExtPtr_ = new labelList(label(sum(scalar(1) - gammaI)), -1);
    labelList& deadCellsExt = *deadCellsExtPtr_;

    label counter = 0;
    forAll (gammaI, cellI)
    {
        if (gammaI[cellI] < SMALL)
        {
            deadCellsExt[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeDeadFaces() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeDeadFaces() const")
            << "create list of dead cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (deadFacesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeDeadFaces() const")
            << "list of dead cells already exist"
            << abort(FatalError);
    }

    deadFacesPtr_ = new labelList(mesh_.nFaces());
    labelList& df = *deadFacesPtr_;
    label nDf = 0;

    const volScalarField& gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    const volScalarField::GeometricBoundaryField& gEPatches =
        gE.boundaryField();

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    forAll (neighbour, faceI)
    {
        if (gammaExtI[neighbour[faceI]] + gammaExtI[owner[faceI]] < SMALL)
        {
            df[nDf] = faceI;
            nDf++;
        }
    }

    forAll (gEPatches, patchI)
    {
        const label start = mesh_.boundaryMesh()[patchI].start();

        if (gEPatches[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gEPatches[patchI].patchInternalField();

            scalarField gammaExtNei =
                gEPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if
                (
                    gammaExtNei[patchFaceI] + gammaExtOwn[patchFaceI] < SMALL
                )
                {
                    df[nDf] = start + patchFaceI;
                    nDf++;
                }
            }
        }
        else
        {
            scalarField gammaExtOwn =
                gEPatches[patchI].patchInternalField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if (gammaExtOwn[patchFaceI] < SMALL)
                {
                    df[nDf] = start + patchFaceI;
                    nDf++;
                }
            }
        }
    }

    df.setSize(nDf);
}


void Foam::immersedBoundaryFvPatch::makeLiveCells() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeLiveCells() const")
            << "create list of live cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (liveCellsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeLiveCells() const")
            << "list of live cells already exist"
            << abort(FatalError);
    }

    const scalarField& gammaI = gamma().internalField();

    liveCellsPtr_ = new labelList(label(sum(gammaI)), -1);
    labelList& liveCells = *liveCellsPtr_;

    label counter = 0;
    forAll (gammaI, cellI)
    {
        if (gammaI[cellI] > (1.0 - SMALL))
        {
            liveCells[counter++] = cellI;
        }
    }
}


void Foam::immersedBoundaryFvPatch::makeIbCellSizes() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbCellsSize() const")
            << "create average sizes of immersed boundary cells"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibCellSizesPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbCellsSize() const")
            << "average sizes of immersed boundary cells already exist"
            << abort(FatalError);
    }

    ibCellSizesPtr_ = new scalarField(ibPoints().size(), 0.0);
    scalarField& delta = *ibCellSizesPtr_;

    if (mesh_.nGeometricD() == 3)
    {
        // Create a list of volumes with mapping to contain only IB cells
        scalarField V(mesh_.V().field(), ibCells());

        delta = Foam::pow(V, 1.0/3.0);
    }
    else
    {
        // For 2-D simulations with the immersed boundary method the geometry
        // needs to be aligned with the z-direction.
        // Having the x- or y-direction as empty is not allowed because of
        // the way the polynomials are expanded
        const Vector<label>& directions = mesh_.geometricD();

        if (directions[0] == -1 || directions[1] == -1)
        {
            FatalErrorIn("immersedBoundaryFvPatch::makeIbCellsSize() const")
                << "For 2-D simulations with the immersed boundary method "
                << "the geometry needs to be aligned with the z-direction.  "
                << "Having the x- or y-direction as empty is not allowed "
                << "because of the way the polynomials are expanded."
                << abort(FatalError);
        }

        scalar thickness = 0.0;

        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        // Field created with mapping for IB cells only
        delta = sqrt(scalarField(mesh_.V().field(), ibCells())/thickness);
    }
}


void Foam::immersedBoundaryFvPatch::makeIbSf() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbSf() const")
            << "creating ibSf and ibMagSf field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibSfPtr_ || ibMagSfPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbSf() const")
            << "ibSf and ibMagSf field already exist"
            << abort(FatalError);
    }

    const vectorField& areas = mesh_.faceAreas();

    // Field created with mapping for IB cells only
    ibSfPtr_ = new vectorField(areas, ibFaces());
    ibMagSfPtr_ = new scalarField(mag(*ibSfPtr_));
}


void Foam::immersedBoundaryFvPatch::makeIbDelta() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeIbDelta() const")
            << "creating delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibDeltaPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeIbDelta() const")
            << "delta field already exist"
            << abort(FatalError);
    }

    const vectorField& C = mesh_.cellCentres();

    // Field created with mapping for IB cells only
    ibDeltaPtr_ =
        new scalarField(mag(ibPoints() - vectorField(C, ibCells())) + SMALL);
}


void Foam::immersedBoundaryFvPatch::makeIbSamplingPointDelta() const
{
    if (debug)
    {
        InfoIn
        (
            "void immersedBoundaryFvPatch::makeIbSamplingPointDelta() const"
        )   << "creating sampling point delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (ibSamplingPointDeltaPtr_)
    {
        FatalErrorIn
        (
            "immersedBoundaryFvPatch::makeIbSamplingPointDelta() const"
        )   << "sampling point delta field already exist"
            << abort(FatalError);
    }

    // Field created with mapping for IB cells only
    ibSamplingPointDeltaPtr_ =
        new scalarField(mag(ibPoints() - ibSamplingPoints()) + SMALL);
}


void Foam::immersedBoundaryFvPatch::makeTriSf() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriSf() const")
            << "creating delta field"
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (triSfPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriSf() const")
            << "triSf field already exist"
            << abort(FatalError);
    }


    const triSurface& triMesh = ibPolyPatch_.ibMesh();
    const pointField& triMeshPoints = triMesh.points();

    triSfPtr_ = new vectorField(triMesh.size());
    vectorField& Sf = *triSfPtr_;

    forAll (triMesh, faceI)
    {
        Sf[faceI] = triMesh[faceI].normal(triMeshPoints);
    }

    if (!ibPolyPatch_.internalFlow())
    {
        // Tri surface points the wrong way; flip all area vectors
        Sf *= -1;
    }
}


void Foam::immersedBoundaryFvPatch::calcCellSearch() const
{
    if (cellSearchPtr_)
    {
        FatalErrorIn("void immersedBoundaryFvPatch::calcCellSearch() const")
            << "Cell tree already calculated"
            << abort(FatalError);
    }

    // Create the octree search for this mesh.  It will be used by other
    // processors when searching for  cells

    // Bounding box containing only local region cells
    treeBoundBox overallBb(mesh_.allPoints());
    Random rndGen(123456);
    overallBb = overallBb.extend(rndGen, 1E-4);
    overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
    overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

    // Collect eligible cells
    const scalarField& gammaExtI = gammaExt().internalField();

    labelList eligibleCells(gammaExtI.size());
    label nElig = 0;

    forAll (gammaExtI, cellI)
    {
        // Note: checking 0:1 indicator
        if (gammaExtI[cellI] > SMALL)
        {
            eligibleCells[nElig] = cellI;
            nElig++;
        }
    }

    // Resize eligible list
    eligibleCells.setSize(nElig);

    if (nElig == 0)
    {
        WarningIn("void immersedBoundaryFvPatch::calcCellSearch() const")
            << "No eligible cells: cellSearch empty"
            << endl;
    }

    // Search
    cellSearchPtr_ = new indexedOctree<treeDataCell>
    (
        treeDataCell
        (
            false,  //  Cache bb.  Reconsider for moving mesh cases
            mesh_,
            eligibleCells
        ),
        overallBb,  // overall search domain
        8,          // maxLevel
        10,         // leafsize
        3.0         // duplicity
    );
}


// Find the cell with the nearest cell centre
void Foam::immersedBoundaryFvPatch::findCellCells
(
    const vector& pt,
    const label cellID,
    labelList& cellCells
) const
{
    const labelListList& cellPoints = mesh_.cellPoints();
    const labelListList& pointCells = mesh_.pointCells();

    const scalarField& gammaExtI = gammaExt().internalField();

    // Create a hash set with estimated size
    labelHashSet cellSet(maxCellCellRows_()*primitiveMesh::facesPerCell_);
    cellSet.insert(cellID);

    // First row
    const labelList& curCellPoints = cellPoints[cellID];

    forAll (curCellPoints, pointI)
    {
        const label curPoint = curCellPoints[pointI];
        const labelList& curPointCells = pointCells[curPoint];

        forAll (curPointCells, cI)
        {
            // Note: checking 0:1 indicator
            if (gammaExtI[curPointCells[cI]] > SMALL)
            {
                // Insert cell as support
                // Search not needed: duplicates automatically filtered
                // HJ, 2/May/2017
                cellSet.insert(curPointCells[cI]);
            }
        }
    }

    labelList curCells = cellSet.toc();

    // Second and other rows
    for (label nRows = 1; nRows < maxCellCellRows_(); nRows++)
    {
        curCells = cellSet.toc();

        forAll (curCells, cellI)
        {
            const label curCell = curCells[cellI];
            const labelList& curCellPoints = cellPoints[curCell];

            forAll (curCellPoints, pointI)
            {
                const label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = pointCells[curPoint];

                forAll (curPointCells, cI)
                {
                    if (gammaExtI[curPointCells[cI]] > SMALL)
                    {
                        // Insert cell as support
                        // Search not needed: duplicates automatically filtered
                        // HJ, 2/May/2017
                        cellSet.insert(curPointCells[cI]);
                    }
                }
            }
        }
    }

    // Erase current cell
    cellSet.erase(cellID);

    // Sorting cells
    const vectorField& C = mesh_.cellCentres();

    curCells = cellSet.toc();
    scalarField distances(curCells.size(), 0);

    forAll (distances, cI)
    {
        distances[cI] =
            mag(C[curCells[cI]] - pt);
    }

    SortableList<scalar> sortedDistances(distances);

    labelList sortedCells(curCells.size(), -1);

    for (label i = 0; i < sortedCells.size(); i++)
    {
        sortedCells[i] =
            curCells[sortedDistances.indices()[i]];
    }

    cellCells = sortedCells;
}


Foam::scalar Foam::immersedBoundaryFvPatch::cellSize(label cellID) const
{
    // Inefficient: reconsider.  HJ, 4/Oct/2017
    if (mesh_.nGeometricD() == 3)
    {
        return Foam::pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0;

        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        return sqrt(mesh_.V().field()[cellID]/thickness);
    }
}


Foam::scalar Foam::immersedBoundaryFvPatch::cellProjection
(
    label cellID,
    const vector& dir
) const
{
    const vectorField& C = mesh_.cellCentres();
    const vectorField& Cf = mesh_.faceCentres();
    const vectorField& Sf = mesh_.faceAreas();

    const labelList& cellFaces = mesh_.cells()[cellID];

    scalar area = 0;

    forAll (cellFaces, faceI)
    {
        label curFace = cellFaces[faceI];

        vector curSf = Sf[curFace];

        if ((curSf & (Cf[curFace] - C[cellID])) < 0)
        {
            curSf *= -1;
        }

        if ((curSf&dir) > 1)
        {
            area += (curSf & dir);
        }
    }

    return area;
}


void Foam::immersedBoundaryFvPatch::clearOut()
{
    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(gammaExtPtr_);
    deleteDemandDrivenData(sGammaPtr_);
    deleteDemandDrivenData(ibCellsPtr_);
    deleteDemandDrivenData(ibFacesPtr_);
    deleteDemandDrivenData(ibFaceCellsPtr_);
    deleteDemandDrivenData(ibFaceFlipsPtr_);
    deleteDemandDrivenData(ibInsideFacesPtr_);
    deleteDemandDrivenData(ibInternalFacesPtr_);
    deleteDemandDrivenData(ibPointsPtr_);
    deleteDemandDrivenData(ibNormalsPtr_);
    deleteDemandDrivenData(hitFacesPtr_);
    deleteDemandDrivenData(ibSamplingPointsPtr_);

    deleteDemandDrivenData(ibSamplingWeightsPtr_);
    deleteDemandDrivenData(ibSamplingProcWeightsPtr_);

    deleteDemandDrivenData(cellsToTriAddrPtr_);
    deleteDemandDrivenData(cellsToTriWeightsPtr_);

    deleteDemandDrivenData(ibCellCellsPtr_);
    deleteDemandDrivenData(ibProcCellsPtr_);
    deleteDemandDrivenData(ibMapPtr_);
    deleteDemandDrivenData(ibProcCentresPtr_);
    deleteDemandDrivenData(ibProcGammaPtr_);
    deleteDemandDrivenData(ibCellProcCellsPtr_);

    deleteDemandDrivenData(deadCellsPtr_);
    deleteDemandDrivenData(deadCellsExtPtr_);
    deleteDemandDrivenData(deadFacesPtr_);
    deleteDemandDrivenData(liveCellsPtr_);
    deleteDemandDrivenData(ibCellSizesPtr_);

    deleteDemandDrivenData(invDirichletMatricesPtr_);
    deleteDemandDrivenData(invNeumannMatricesPtr_);

    deleteDemandDrivenData(ibSfPtr_);
    deleteDemandDrivenData(ibMagSfPtr_);
    deleteDemandDrivenData(ibDeltaPtr_);
    deleteDemandDrivenData(ibSamplingPointDeltaPtr_);

    deleteDemandDrivenData(triSfPtr_);

    deleteDemandDrivenData(cellSearchPtr_);
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::initMovePoints()
{}


void Foam::immersedBoundaryFvPatch::movePoints()
{
    clearOut();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::immersedBoundaryFvPatch::immersedBoundaryFvPatch
(
    const polyPatch& patch,
    const fvBoundaryMesh& bm
)
:
    fvPatch(patch, bm),
    ibPolyPatch_(refCast<const immersedBoundaryPolyPatch>(patch)),
    mesh_(bm.mesh()),
    gammaPtr_(NULL),
    gammaExtPtr_(NULL),
    sGammaPtr_(NULL),
    ibCellsPtr_(NULL),
    ibFacesPtr_(NULL),
    ibFaceCellsPtr_(NULL),
    ibFaceFlipsPtr_(NULL),
    ibInsideFacesPtr_(NULL),
    ibInternalFacesPtr_(NULL),
    ibPointsPtr_(NULL),
    ibNormalsPtr_(NULL),
    hitFacesPtr_(NULL),
    ibSamplingPointsPtr_(NULL),
    ibSamplingWeightsPtr_(NULL),
    ibSamplingProcWeightsPtr_(NULL),
    cellsToTriAddrPtr_(NULL),
    cellsToTriWeightsPtr_(NULL),
    ibCellCellsPtr_(NULL),
    ibProcCellsPtr_(NULL),
    ibMapPtr_(NULL),
    ibProcCentresPtr_(NULL),
    ibProcGammaPtr_(NULL),
    ibCellProcCellsPtr_(NULL),
    deadCellsPtr_(NULL),
    deadCellsExtPtr_(NULL),
    deadFacesPtr_(NULL),
    liveCellsPtr_(NULL),
    ibCellSizesPtr_(NULL),
    invDirichletMatricesPtr_(NULL),
    invNeumannMatricesPtr_(NULL),
    ibSfPtr_(NULL),
    ibMagSfPtr_(NULL),
    ibDeltaPtr_(NULL),
    ibSamplingPointDeltaPtr_(NULL),
    triSfPtr_(NULL),
    cellSearchPtr_(NULL)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::volScalarField& Foam::immersedBoundaryFvPatch::gamma() const
{
    if (!gammaPtr_)
    {
        makeGamma();
    }

    return *gammaPtr_;
}


const Foam::volScalarField& Foam::immersedBoundaryFvPatch::gammaExt() const
{
    if (!gammaExtPtr_)
    {
        makeGammaExt();
    }

    return *gammaExtPtr_;
}


const Foam::surfaceScalarField& Foam::immersedBoundaryFvPatch::sGamma() const
{
    if (!sGammaPtr_)
    {
        makeSGamma();
    }

    return *sGammaPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibCells() const
{
    if (!ibCellsPtr_)
    {
        makeIbCells();
    }

    return *ibCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibFaces() const
{
    if (!ibFacesPtr_)
    {
        makeIbFaces();
    }

    return *ibFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibFaceCells() const
{
    if (!ibFaceCellsPtr_)
    {
        makeIbFaces();
    }

    return *ibFaceCellsPtr_;
}


const Foam::boolList& Foam::immersedBoundaryFvPatch::ibFaceFlips() const
{
    if (!ibFaceFlipsPtr_)
    {
        makeIbFaces();
    }

    return *ibFaceFlipsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibInsideFaces() const
{
    if (!ibInsideFacesPtr_)
    {
        makeIbInsideFaces();
    }

    return *ibInsideFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::ibInternalFaces() const
{
    if (!ibInternalFacesPtr_)
    {
        makeIbInternalFaces();
    }

    return *ibInternalFacesPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibPoints() const
{
    if (!ibPointsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibPointsPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibNormals() const
{
    if (!ibNormalsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibNormalsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::hitFaces() const
{
    if (!hitFacesPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *hitFacesPtr_;
}


const Foam::vectorField&
Foam::immersedBoundaryFvPatch::ibSamplingPoints() const
{
    if (!ibSamplingPointsPtr_)
    {
        makeIbPointsAndNormals();
    }

    return *ibSamplingPointsPtr_;
}


const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibCellCells() const
{
    if (!ibCellCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibCellCellsPtr_;
}


const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibProcCells() const
{
    if (!ibProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCellsPtr_;
}


const Foam::mapDistribute&
Foam::immersedBoundaryFvPatch::ibMap() const
{
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const mapDistribute& immersedBoundaryFvPatch::ibMap() const"
        )   << "Requested mapDistribute in serial run.  "
            << "This is not allowed"
            << abort(FatalError);
    }

    if (!ibCellProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibMapPtr_;
}


const Foam::vectorList&
Foam::immersedBoundaryFvPatch::ibProcCentres() const
{
    if (!ibProcCentresPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCentresPtr_;
}


const Foam::scalarList&
Foam::immersedBoundaryFvPatch::ibProcGamma() const
{
    if (!ibProcGammaPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcGammaPtr_;
}


const Foam::labelListList&
Foam::immersedBoundaryFvPatch::ibCellProcCells() const
{
    if (!ibCellProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibCellProcCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadCells() const
{
    if (!deadCellsPtr_)
    {
        makeDeadCells();
    }

    return *deadCellsPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadCellsExt() const
{
    if (!deadCellsExtPtr_)
    {
        makeDeadCellsExt();
    }

    return *deadCellsExtPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::deadFaces() const
{
    if (!deadFacesPtr_)
    {
        makeDeadFaces();
    }

    return *deadFacesPtr_;
}


const Foam::labelList& Foam::immersedBoundaryFvPatch::liveCells() const
{
    if (!liveCellsPtr_)
    {
        makeLiveCells();
    }

    return *liveCellsPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibCellSizes() const
{
    if (!ibCellSizesPtr_)
    {
        makeIbCellSizes();
    }

    return *ibCellSizesPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::ibSf() const
{
    if (!ibSfPtr_)
    {
        makeIbSf();
    }

    return *ibSfPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibMagSf() const
{
    if (!ibMagSfPtr_)
    {
        makeIbSf();
    }

    return *ibMagSfPtr_;
}


const Foam::scalarField& Foam::immersedBoundaryFvPatch::ibDelta() const
{
    if (!ibDeltaPtr_)
    {
        makeIbDelta();
    }

    return *ibDeltaPtr_;
}


const Foam::scalarField&
Foam::immersedBoundaryFvPatch::ibSamplingPointDelta() const
{
    if (!ibSamplingPointDeltaPtr_)
    {
        makeIbSamplingPointDelta();
    }

    return *ibSamplingPointDeltaPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::triSf() const
{
    if (!triSfPtr_)
    {
        makeTriSf();
    }

    return *triSfPtr_;
}


const Foam::vectorField& Foam::immersedBoundaryFvPatch::triCf() const
{
    return ibMesh().faceCentres();
}


const Foam::indexedOctree<Foam::treeDataCell>&
Foam::immersedBoundaryFvPatch::cellSearch() const
{
    if (!cellSearchPtr_)
    {
        calcCellSearch();
    }

    return *cellSearchPtr_;
}


const Foam::dynamicLabelList&
Foam::immersedBoundaryFvPatch::triFacesInMesh() const
{
    // Check if triFacesInMesh has been updated this time step
    if (ibUpdateTimeIndex_ != mesh_.time().timeIndex())
    {
        const vectorField& triCf = this->triCf();

        triFacesInMesh_.clear();
        triFacesInMesh_.setCapacity(triCf.size()/2);

        // Use octree search to find the cell
        treeBoundBox overallBb(mesh_.points());
        Random rndGen(123456);
        overallBb = overallBb.extend(rndGen, 1E-4);
        overallBb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
        overallBb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

        // Search
        indexedOctree<treeDataCell> cellSearch
        (
            treeDataCell(false, mesh_),
            overallBb,  // overall search domain
            8,          // maxLevel
            10,         // leafsize
            3.0         // duplicity
        );
        scalar span = cellSearch.bb().mag();

        // Find tri faces with centre inside the processor mesh
        forAll (triCf, fI)
        {
            const vector curTriCf = triCf[fI];

            if (!mesh_.bounds().containsInside(curTriCf))
            {
                // Face centre is not inside the mesh
                continue;
            }
            else
            {
#if 0
                // Slow mesh cell search
                if (mesh_.findCell(curTriCf) != -1)
                {
                    triFacesInMesh_.append(fI);
                }
#else
                // Octree mesh search
                pointIndexHit pih = cellSearch.findNearest(curTriCf, span);

                if (pih.hit())
                {
                    // Found a hit.  Additional check for point in cell
                    const label hitCell = pih.index();

                    if
                    (
                        mesh_.pointInCellBB
                        (
                            curTriCf,
                            hitCell
                        )
                    )
                    {
                        triFacesInMesh_.append(fI);
                    }
                }
#endif
            }
        }

        ibUpdateTimeIndex_ = mesh_.time().timeIndex();
    }

    return triFacesInMesh_;
}


// ************************************************************************* //
