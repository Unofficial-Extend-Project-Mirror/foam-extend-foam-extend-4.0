/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
//     5
);


const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::angleFactor_
(
    "immersedBoundaryAngleFactor",
    80
//     170
);


const Foam::debug::optimisationSwitch
Foam::immersedBoundaryFvPatch::maxCellCellRows_
(
    "immersedBoundaryMaxCellCellRows",
    4
//     5
);


const Foam::debug::tolerancesSwitch
Foam::immersedBoundaryFvPatch::distFactor_
(
    "immersedBoundaryDistFactor",
    1.5
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

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

    forAll (ibc, cellI)
    {
        gammaI[ibc[cellI]] = 0.0;
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

    labelHashSet ibCellSet;

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    forAll (neighbour, faceI)
    {
        if (mag(gammaExtI[neighbour[faceI]] - gammaExtI[owner[faceI]]) > SMALL)
        {
            if (gammaExtI[owner[faceI]] > SMALL)
            {
                if (!ibCellSet.found(owner[faceI]))
                {
                    ibCellSet.insert(owner[faceI]);
                }
            }
            else
            {
                if (!ibCellSet.found(neighbour[faceI]))
                {
                    ibCellSet.insert(neighbour[faceI]);
                }
            }
        }
    }

    forAll (gE.boundaryField(), patchI)
    {
        if (gE.boundaryField()[patchI].coupled())
        {
            scalarField gammaExtOwn =
                gE.boundaryField()[patchI].patchInternalField();

            scalarField gammaExtNei =
                gE.boundaryField()[patchI].patchNeighbourField();

            const unallocLabelList& fCells =
                mesh_.boundary()[patchI].faceCells();

            forAll (gammaExtOwn, faceI)
            {
                if
                (
                    mag(gammaExtNei[faceI] - gammaExtOwn[faceI])
                  > SMALL
                )
                {
                    if (gammaExtOwn[faceI] > SMALL)
                    {
                        if (!ibCellSet.found(fCells[faceI]))
                        {
                            ibCellSet.insert(fCells[faceI]);
                        }
                    }
                    else if (2*gammaExtOwn.size() == fCells.size())
                    {
                        if
                        (
                           !ibCellSet.found
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            )
                        )
                        {
                            ibCellSet.insert
                            (
                                fCells[gammaExtOwn.size() + faceI]
                            );
                        }
                    }
                }
            }
        }
    }

    ibCellsPtr_ = new labelList(ibCellSet.toc());
    sort(*ibCellsPtr_);

    Pout << "Number of IB cells: " << ibCellsPtr_->size() << endl;
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
    labelList cornerCells;

    const triSurfaceSearch& tss = ibPolyPatch_.triSurfSearch();

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

        labelHashSet cornerIbCellSet;

        const labelList& ibf = ibFaces();

        forAll (ibf, faceI)
        {
            const label& ownCell = own[ibf[faceI]];
            const label& neiCell = nei[ibf[faceI]];

//             label ibCell = -1;
            label liveCell = -1;

            if (gammaI[ownCell] > SMALL)
            {
//                 ibCell = neiCell;
                liveCell = ownCell;
            }
            else
            {
//                 ibCell = ownCell;
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

                        //HJ bug fix?
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
                    if (!cornerIbCellSet.found(liveCell))
                    {
                        cornerIbCellSet.insert(liveCell);
                    }
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

    DynamicList<label> ibF(2*ibc.size());
    DynamicList<label> ibFC(2*ibc.size());
    DynamicList<bool> ibFF(2*ibc.size());

    const unallocLabelList& owner = mesh_.owner();
    const unallocLabelList& neighbour = mesh_.neighbour();

//     const scalarField& gammaI = gamma().internalField();

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

    labelHashSet ibInsideFSet;

    const labelList& owner = mesh_.faceOwner();
    const labelList& neighbour = mesh_.faceNeighbour();

    const volScalarField gE = gammaExt();
    const scalarField& gammaExtI = gE.internalField();

    forAll (neighbour, faceI)
    {
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

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    mag(gammaNei[faceI] - gammaOwn[faceI]) > SMALL
                )
                {
                    if (!ibInsideFSet.found(start + faceI))
                    {
                        ibInsideFSet.insert(start + faceI);
                    }

                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInsideFSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInsideFSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    ibInsideFacesPtr_ = new labelList(ibInsideFSet.toc());
    sort(*ibInsideFacesPtr_);
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

    labelHashSet ibInternalFacesSet;

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

            label size = mesh_.boundaryMesh()[patchI].size();
            label start = mesh_.boundaryMesh()[patchI].start();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    (gammaNei[faceI] > SMALL)
                 && (gammaOwn[faceI] > SMALL)
                )
                {
                    if (!ibInternalFacesSet.found(start + faceI))
                    {
                        ibInternalFacesSet.insert(start + faceI);
                    }

                    if (2*gammaOwn.size() == size)
                    {
                        if
                        (
                           !ibInternalFacesSet.found
                            (
                                start + size/2 + faceI
                            )
                        )
                        {
                            ibInternalFacesSet.insert
                            (
                                start + size/2 + faceI
                            );
                        }
                    }
                }
            }
        }
    }

    ibInternalFacesPtr_ = new labelList(ibInternalFacesSet.toc());
    sort(*ibInternalFacesPtr_);
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

    // Find average cell dimension
    const labelList& ibc = ibCells();

    scalarField delta(ibc.size(), 0.0);

    forAll (delta, cellI)
    {
        delta[cellI] = cellSize(ibc[cellI]);
    }

    // Find nearest triSurface point for each interface cell centre
    ibPointsPtr_ = new vectorField(ibc.size(), vector::zero);
    ibNormalsPtr_ = new vectorField(ibc.size(), vector::zero);
    hitFacesPtr_ = new labelList(ibc.size(), -1);
    ibSamplingPointsPtr_ = new vectorField(ibc.size(), vector::zero);

    vectorField& ibPoints = *ibPointsPtr_;
    vectorField& ibNormals = *ibNormalsPtr_;
    labelList& ibHitFaces = *hitFacesPtr_;
    vectorField& ibSamplingPoints = *ibSamplingPointsPtr_;

    const vectorField& C = mesh_.C().internalField();

    // Get IB cell centres
    vectorField ibCellCentres(C, ibc);

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
                << "Hit data = " << pih << nl
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

    const labelList& ibc = ibCells();

    ibCellCellsPtr_ = new labelListList(ibc.size());
    labelListList& cellCells = *ibCellCellsPtr_;

    ibProcCentresPtr_ = new FieldField<Field, vector>(Pstream::nProcs());
    FieldField<Field, vector>& procCentres = *ibProcCentresPtr_;

    forAll (procCentres, procI)
    {
        procCentres.set(procI, new vectorField(0));
    }

    ibProcGammaPtr_ = new FieldField<Field, scalar>(Pstream::nProcs());
    FieldField<Field, scalar>& procGamma = *ibProcGammaPtr_;

    forAll (procGamma, procI)
    {
        procGamma.set(procI, new scalarField(0));
    }


    ibCellProcCellsPtr_ = new List<List<labelPair> >(ibc.size());
    List<List<labelPair> >& cellProcCells = *ibCellProcCellsPtr_;

    const cellList& meshCells = mesh_.cells();

    // In geometry initialisation, fields are not available: use raw mesh data
    // HJ after ZT, 6/Dec/2012
    const vectorField& C = mesh_.cellCentres();

    scalarField rM(ibCellSizes());
    rM *= radiusFactor_();

    const vectorField& ibp = ibPoints();

    // Note: the algorithm is originally written with inward-facing normals
    // and subsequently changed: IB surface normals point outwards
    // HJ, 21/May/2012
//     const vectorField& ibn = ibNormals();

    forAll (cellCells, cellI)
    {
        labelList curCells;

        findCellCells
        (
            C[ibc[cellI]],
            ibc[cellI],
            curCells
        );

        cellCells[cellI] = labelList(curCells.size(), -1);

        label cI = 0;

        for (label i = 0; i < curCells.size(); i++)
        {
            label curCell = curCells[i];
            scalar r = mag(C[curCell] - C[ibc[cellI]]);

            // Collect the cells within rM of the fitting cell
            if (r <= rM[cellI])
            {
//                 scalar angleLimit =
//                     -Foam::cos(angleFactor_()*mathematicalConstant::pi/180);

                vector dir = (C[curCell] - ibp[cellI]);
                dir /= mag(dir) + SMALL;

                // Change of sign of normal.  HJ, 21/May/2012
//                 if ((-ibn[cellI] & dir) >= angleLimit)
                {
                    cellCells[cellI][cI++] = curCell;
                }
            }
        }

        cellCells[cellI].setSize(cI);
    }

    if (Pstream::parRun())
    {
        // Find immersed boundary cells next to processor boundaries
        labelHashSet procIbCellsSet;

        forAll (ibc, cellI)
        {
            const labelList& curCellCells = cellCells[cellI];

            if (curCellCells.size())
            {
                forAll (curCellCells, cI)
                {
                    const labelList& faces = meshCells[curCellCells[cI]];

                    bool foundProcessorFace = false;

                    forAll (faces, faceI)
                    {
                        label patchID =
                            mesh_.boundaryMesh().whichPatch(faces[faceI]);

                        if (patchID != -1)
                        {
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
                const labelList& faces = meshCells[ibc[cellI]];

                bool foundProcessorFace = false;

                forAll (faces, faceI)
                {
                    label patchID =
                        mesh_.boundaryMesh().whichPatch(faces[faceI]);

                    if (patchID != -1)
                    {
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
        labelList procIbCells = procIbCellsSet.toc();
        sort(procIbCells);

        // Note: consider more sophisticated gather-scatter
        // HJ, 18/Jun/2015

        // Send and receive number  of immersed boundary cells
        // next to processor boundaries
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << procIbCells.size();
                }
            }
        }

        labelList sizes(Pstream::nProcs(), 0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> sizes[procI];
                }
            }
        }

        // Send and receive ibc centres and radii
        vectorField centres(procIbCells.size(), vector::zero);
        forAll (centres, cellI)
        {
            centres[cellI] = C[ibc[procIbCells[cellI]]];
        }
        scalarField procRMax(rM, procIbCells);

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        centres.size()*sizeof(vector)
                      + procRMax.size()*sizeof(scalar)
                    );

                    toProc << centres << procRMax;
                }
            }
        }

        FieldField<Field, vector> ctrs(Pstream::nProcs());
        FieldField<Field, scalar> rMax(Pstream::nProcs());

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                ctrs.set
                (
                    procI,
                    new vectorField(sizes[procI], vector::zero)
                );

                rMax.set
                (
                    procI,
                    new scalarField(sizes[procI], 0)
                );

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        sizes[procI]*sizeof(vector)
                      + sizes[procI]*sizeof(scalar)
                    );

                    fromProc >> ctrs[procI] >> rMax[procI];
                }
            }
            else
            {
                ctrs.set(procI, new vectorField(0));
                rMax.set(procI, new scalarField(0));
            }
        }

        // Find cells needed by other processors
        if (ibProcCellsPtr_)
        {
            FatalErrorIn("immersedBoundaryFvPatch::makeIbCellCells() const")
                << "procCells addressing already exists"
                << abort(FatalError);
        }

        ibProcCellsPtr_ = new labelListList(Pstream::nProcs());
        labelListList& procCells = *ibProcCellsPtr_;

        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            labelHashSet procCellSet;

            if (procI != Pstream::myProcNo())
            {
                forAll (ctrs[procI], cellI)
                {
                    label nearestCellID =
                        findNearestCell(ctrs[procI][cellI]);

                    if (nearestCellID == -1)
                    {
                        FatalErrorIn
                        (
                            "immersedBoundaryFvPatch::makeIbCellCells() const"
                        ) << "Can't find nearest cell."
                            << abort(FatalError);
                    }

                    scalar R = mag(C[nearestCellID] - ctrs[procI][cellI]);

                    if (R < rMax[procI][cellI])
                    {
                        if (!procCellSet.found(nearestCellID))
                        {
                            procCellSet.insert(nearestCellID);
                        }

                        labelList tmpCellList;

                        findCellCells
                        (
                            ctrs[procI][cellI],
                            nearestCellID,
                            tmpCellList
                        );

                        forAll (tmpCellList, cI)
                        {
                            scalar r =
                                mag
                                (
                                    C[tmpCellList[cI]]
                                  - ctrs[procI][cellI]
                                );

                            if (r <= rMax[procI][cellI])
                            {
                                if (!procCellSet.found(tmpCellList[cI]))
                                {
                                    procCellSet.insert(tmpCellList[cI]);
                                }
                            }
                        }
                    }
                }
            }

            procCells[procI] = procCellSet.toc();
        }


        // Send and receive sizes
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    toProc << procCells[procI].size();
                }
            }
        }

        labelList procSizes(Pstream::nProcs(), 0);
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        sizeof(label)
                    );

                    fromProc >> procSizes[procI];
                }
            }
        }

        // Send cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                vectorField centres(C, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        centres.size()*sizeof(vector)
                    );

                    toProc << centres;
                }
            }
        }

        // Receive cell centres
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procCentres[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        procSizes[procI]*sizeof(vector)
                    );

                    fromProc >> procCentres[procI];
                }
            }
            // else: already set to zero-size field
        }

        // Send cell gamma
        const scalarField& gammaI = gamma().internalField();
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                scalarField gamma(gammaI, procCells[procI]);

                // Parallel data exchange
                {
                    OPstream toProc
                    (
                        Pstream::blocking,
                        procI,
                        gamma.size()*sizeof(scalar)
                    );

                    toProc << gamma;
                }
            }
        }

        // Receive cell gamma
        for (label procI = 0; procI < Pstream::nProcs(); procI++)
        {
            if (procI != Pstream::myProcNo())
            {
                procGamma[procI].setSize(procSizes[procI]);

                // Parallel data exchange
                {
                    IPstream fromProc
                    (
                        Pstream::blocking,
                        procI,
                        procSizes[procI]*sizeof(scalar)
                    );

                    fromProc >> procGamma[procI];
                }
            }
            // else: already set to zero-size field
        }

        // Cell-procCells addressing
        forAll (cellProcCells, cellI)
        {
            scalar rMax = rM[cellI];

            cellProcCells[cellI].setSize(100);

            label index = 0;
            forAll (procCentres, procI)
            {
                forAll (procCentres[procI], pointI)
                {
                    scalar r =
                        mag
                        (
                            procCentres[procI][pointI]
                          - C[ibc[cellI]]
                        );

                    if (r <= rMax)
                    {
                        cellProcCells[cellI][index].first() = procI;
                        cellProcCells[cellI][index].second() = pointI;
                        index++;
                    }
                }
            }

            cellProcCells[cellI].setSize(index);
        }
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


// Find the cell with the nearest cell centre
Foam::label Foam::immersedBoundaryFvPatch::findNearestCell
(
    const point& location
) const
{
    const vectorField& C = mesh_.cellCentres();

    const scalarField& gammaExtI = gammaExt().internalField();

    label nearestCellI = -1;
    scalar minProximity = GREAT;

    for (label cellI = 0; cellI < C.size(); cellI++)
    {
        if (gammaExtI[cellI] > SMALL)
        {
            scalar proximity = magSqr(C[cellI] - location);

            if (proximity < minProximity)
            {
                nearestCellI = cellI;
                minProximity = proximity;
            }
        }
    }

    return nearestCellI;
}


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

    labelHashSet cellSet;
    cellSet.insert(cellID);

    // First row
    const labelList& curCellPoints = cellPoints[cellID];

    forAll (curCellPoints, pointI)
    {
        label curPoint = curCellPoints[pointI];
        const labelList& curPointCells = pointCells[curPoint];

        forAll (curPointCells, cI)
        {
            if (gammaExtI[curPointCells[cI]] > SMALL)
            {
                if (!cellSet.found(curPointCells[cI]))
                {
                    cellSet.insert(curPointCells[cI]);
                }
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
            label curCell = curCells[cellI];
            const labelList& curCellPoints = cellPoints[curCell];

            forAll (curCellPoints, pointI)
            {
                label curPoint = curCellPoints[pointI];
                const labelList& curPointCells = pointCells[curPoint];

                forAll (curPointCells, cI)
                {
                    if (gammaExtI[curPointCells[cI]] > SMALL)
                    {
                        if (!cellSet.found(curPointCells[cI]))
                        {
                            cellSet.insert(curPointCells[cI]);
                        }
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
    scalar delta;

    if (mesh_.nGeometricD() == 3)
    {
        delta = Foam::pow(mesh_.V().field()[cellID], 1.0/3.0);
    }
    else
    {
        scalar thickness = 0.0;
        const Vector<label>& directions = mesh_.geometricD();
        for (direction dir = 0; dir < directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        delta = sqrt(mesh_.V().field()[cellID]/thickness);
    }

    return delta;
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
    triSfPtr_(NULL)
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


const Foam::FieldField<Foam::Field, Foam::vector>&
Foam::immersedBoundaryFvPatch::ibProcCentres() const
{
    if (!ibProcCentresPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCentresPtr_;
}


const Foam::FieldField<Foam::Field, Foam::scalar>&
Foam::immersedBoundaryFvPatch::ibProcGamma() const
{
    if (!ibProcGammaPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcGammaPtr_;
}


const Foam::List<Foam::List<Foam::labelPair> >&
Foam::immersedBoundaryFvPatch::ibCellProcCells() const
{
    if (!ibCellProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibCellProcCellsPtr_;
}


const Foam::labelListList& Foam::immersedBoundaryFvPatch::ibProcCells() const
{
    if (!ibProcCellsPtr_)
    {
        makeIbCellCells();
    }

    return *ibProcCellsPtr_;
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


// ************************************************************************* //
