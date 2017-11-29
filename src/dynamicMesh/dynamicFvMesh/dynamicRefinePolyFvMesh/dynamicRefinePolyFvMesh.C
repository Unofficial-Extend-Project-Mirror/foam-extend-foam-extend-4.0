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

#include "dynamicRefinePolyFvMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvCFD.H"
#include "syncTools.H"
#include "pointFields.H"
#include "directTopoChange.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(dynamicRefinePolyFvMesh, 0);

addToRunTimeSelectionTable(dynamicFvMesh, dynamicRefinePolyFvMesh, IOobject);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

label dynamicRefinePolyFvMesh::count
(
    const PackedBoolList& l,
    const unsigned int val
)
{
    label n = 0;
    forAll(l, i)
    {
        if (l.get(i) == val)
        {
            n++;
        }

        // debug also serves to get-around Clang compiler trying to optimise
        // out this forAll loop under O3 optimisation

        if (debug)
        {
            Info<< "n=" << n << endl;
        }
    }

    return n;
}


void dynamicRefinePolyFvMesh::readDict()
{
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    correctFluxes_ = List<Pair<word> >(refineDict.lookup("correctFluxes"));

    dumpLevel_ = Switch(refineDict.lookup("dumpLevel"));
}


// Refines cells, maps fields and recalculates (an approximate) flux
autoPtr<mapPolyMesh> dynamicRefinePolyFvMesh::refine
(
    const labelList& cellsToRefine
)
{
    // Mesh changing engine.
    directTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setRefinement(cellsToRefine, meshMod);

    // Create mesh (with inflation), return map from old to new mesh.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Refined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells." << endl;

    if (debug)
    {
        // Check map.
        for (label faceI = 0; faceI < nInternalFaces(); faceI++)
        {
            label oldFaceI = map().faceMap()[faceI];

            if (oldFaceI >= nInternalFaces())
            {
                FatalErrorIn("dynamicRefinePolyFvMesh::refine(const labelList&)")
                    << "New internal face:" << faceI
                    << " fc:" << faceCentres()[faceI]
                    << " originates from boundary oldFace:" << oldFaceI
                    << abort(FatalError);
            }
        }
    }


    // Update fields
    updateMesh(map);

    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified/added faces. All the faces which only
    // have been renumbered will already have been handled by the mapping.
    {
        const labelList& faceMap = map().faceMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        // Storage for any master faces. These will be the original faces on the
        // coarse cell that gets split into four for hex cell and quad
        // face. Rather, the master face gets modified and three faces get added
        // from the master.
        labelHashSet masterFaces(4*cellsToRefine.size());

        forAll(faceMap, faceI)
        {
            label oldFaceI = faceMap[faceI];

            if (oldFaceI >= 0)
            {
                label masterFaceI = reverseFaceMap[oldFaceI];

                if (masterFaceI < 0)
                {
                    FatalErrorIn
                    (
                        "dynamicRefinePolyFvMesh::refine(const labelList&)"
                    )   << "Problem: should not have removed faces"
                        << " when refining."
                        << nl << "face:" << faceI << abort(FatalError);
                }
                else if (masterFaceI != faceI)
                {
                    masterFaces.insert(masterFaceI);
                }
            }
        }
        if (debug)
        {
            Info<< "Found " << returnReduce(masterFaces.size(), sumOp<label>())
                << " split faces " << endl;
        }

        forAll(correctFluxes_, i)
        {
            if (debug)
            {
                Info<< "Mapping flux " << correctFluxes_[i][0]
                    << " using interpolated flux " << correctFluxes_[i][1]
                    << endl;
            }
            surfaceScalarField& phi = const_cast<surfaceScalarField&>
            (
                lookupObject<surfaceScalarField>(correctFluxes_[i][0])
            );
            surfaceScalarField phiU =
                fvc::interpolate
                (
                    lookupObject<volVectorField>(correctFluxes_[i][1])
                )
              & Sf();

            // Recalculate new internal faces.
            for (label faceI = 0; faceI < nInternalFaces(); faceI++)
            {
                label oldFaceI = faceMap[faceI];

                if (oldFaceI == -1)
                {
                    // Inflated/appended
                    phi[faceI] = phiU[faceI];
                }
                else if (reverseFaceMap[oldFaceI] != faceI)
                {
                    // face-from-masterface
                    phi[faceI] = phiU[faceI];
                }
            }

            // Recalculate new boundary faces.
            forAll(phi.boundaryField(), patchI)
            {
                fvsPatchScalarField& patchPhi = phi.boundaryField()[patchI];
                const fvsPatchScalarField& patchPhiU =
                    phiU.boundaryField()[patchI];

                label faceI = patchPhi.patch().patch().start();

                forAll(patchPhi, i)
                {
                    label oldFaceI = faceMap[faceI];

                    if (oldFaceI == -1)
                    {
                        // Inflated/appended
                        patchPhi[i] = patchPhiU[i];
                    }
                    else if (reverseFaceMap[oldFaceI] != faceI)
                    {
                        // face-from-masterface
                        patchPhi[i] = patchPhiU[i];
                    }

                    faceI++;
                }
            }

            // Update master faces
            forAllConstIter(labelHashSet, masterFaces, iter)
            {
                label faceI = iter.key();

                if (isInternalFace(faceI))
                {
                    phi[faceI] = phiU[faceI];
                }
                else
                {
                    label patchI = boundaryMesh().whichPatch(faceI);
                    label i = faceI - boundaryMesh()[patchI].start();

                    const fvsPatchScalarField& patchPhiU =
                        phiU.boundaryField()[patchI];

                    fvsPatchScalarField& patchPhi =
                        phi.boundaryField()[patchI];

                    patchPhi[i] = patchPhiU[i];
                }
            }
        }
    }



    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


// Combines previously split cells, maps fields and recalculates
// (an approximate) flux
autoPtr<mapPolyMesh> dynamicRefinePolyFvMesh::unrefine
(
    const labelList& splitPoints
)
{
    directTopoChange meshMod(*this);

    // Play refinement commands into mesh changer.
    meshCutter_.setUnrefinement(splitPoints, meshMod);


    // Save information on faces that will be combined
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Find the faceMidPoints on cells to be combined.
    // for each face resulting of split of face into four store the
    // midpoint
    Map<label> faceToSplitPoint(3*splitPoints.size());

    {
        forAll(splitPoints, i)
        {
            label pointI = splitPoints[i];

            const labelList& pEdges = pointEdges()[pointI];

            forAll(pEdges, j)
            {
                label otherPointI = edges()[pEdges[j]].otherVertex(pointI);

                const labelList& pFaces = pointFaces()[otherPointI];

                forAll(pFaces, pFaceI)
                {
                    faceToSplitPoint.insert(pFaces[pFaceI], otherPointI);
                }
            }
        }
    }


    // Change mesh and generate map.
    //autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, true);
    autoPtr<mapPolyMesh> map = meshMod.changeMesh(*this, false);

    Info<< "Unrefined from "
        << returnReduce(map().nOldCells(), sumOp<label>())
        << " to " << globalData().nTotalCells() << " cells."
        << endl;

    // Update fields
    updateMesh(map);


    // Move mesh
    /*
    pointField newPoints;
    if (map().hasMotionPoints())
    {
        newPoints = map().preMotionPoints();
    }
    else
    {
        newPoints = points();
    }
    movePoints(newPoints);
    */

    // Correct the flux for modified faces.
    {
        const labelList& reversePointMap = map().reversePointMap();
        const labelList& reverseFaceMap = map().reverseFaceMap();

        forAll(correctFluxes_, i)
        {
            if (debug)
            {
                Info<< "Mapping flux " << correctFluxes_[i][0]
                    << " using interpolated flux " << correctFluxes_[i][1]
                    << endl;
            }
            surfaceScalarField& phi = const_cast<surfaceScalarField&>
            (
                lookupObject<surfaceScalarField>(correctFluxes_[i][0])
            );
            surfaceScalarField phiU =
                fvc::interpolate
                (
                    lookupObject<volVectorField>(correctFluxes_[i][1])
                )
              & Sf();

            forAllConstIter(Map<label>, faceToSplitPoint, iter)
            {
                label oldFaceI = iter.key();
                label oldPointI = iter();

                if (reversePointMap[oldPointI] < 0)
                {
                    // midpoint was removed. See if face still exists.
                    label faceI = reverseFaceMap[oldFaceI];

                    if (faceI >= 0)
                    {
                        if (isInternalFace(faceI))
                        {
                            phi[faceI] = phiU[faceI];
                        }
                        else
                        {
                            label patchI = boundaryMesh().whichPatch(faceI);
                            label i = faceI - boundaryMesh()[patchI].start();

                            const fvsPatchScalarField& patchPhiU =
                                phiU.boundaryField()[patchI];

                            fvsPatchScalarField& patchPhi =
                                phi.boundaryField()[patchI];

                            patchPhi[i] = patchPhiU[i];
                        }
                    }
                }
            }
        }
    }


    // Update numbering of cells/vertices.
    meshCutter_.updateMesh(map);

    // Debug: Check refinement levels (across faces only)
    meshCutter_.checkRefinementLevels(-1, labelList(0));

    return map;
}


// Get max of connected point
scalarField dynamicRefinePolyFvMesh::maxPointField(const scalarField& pFld) const
{
    scalarField vFld(nCells(), -GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            vFld[pCells[i]] = max(vFld[pCells[i]], pFld[pointI]);
        }
    }
    return vFld;
}


// Get min of connected cell
scalarField dynamicRefinePolyFvMesh::minCellField(const volScalarField& vFld) const
{
    scalarField pFld(nPoints(), GREAT);

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        forAll(pCells, i)
        {
            pFld[pointI] = min(pFld[pointI], vFld[pCells[i]]);
        }
    }
    return pFld;
}


// Simple (non-parallel) interpolation by averaging.
scalarField dynamicRefinePolyFvMesh::cellToPoint(const scalarField& vFld) const
{
    scalarField pFld(nPoints());

    forAll(pointCells(), pointI)
    {
        const labelList& pCells = pointCells()[pointI];

        scalar sum = 0.0;
        forAll(pCells, i)
        {
            sum += vFld[pCells[i]];
        }
        pFld[pointI] = sum/pCells.size();
    }
    return pFld;
}


// Calculate error. Is < 0 or distance from inbetween levels
scalarField dynamicRefinePolyFvMesh::error
(
    const scalarField& fld,
    const scalar minLevel,
    const scalar maxLevel
) const
{
    const scalar halfLevel = 0.5*(minLevel + maxLevel);

    scalarField c(fld.size(), -1);

    forAll(fld, i)
    {
        if (fld[i] >= minLevel && fld[i] < maxLevel)
        {
            c[i] = mag(fld[i] - halfLevel);
        }
    }
    return c;
}


void dynamicRefinePolyFvMesh::selectRefineCandidates
(
    const scalar lowerRefineLevel,
    const scalar upperRefineLevel,
    const scalarField& vFld,
    PackedBoolList& candidateCell
) const
{
    // Get error per cell. Is -1 (not to be refined) to >0 (to be refined,
    // higher more desirable to be refined).
    scalarField cellError
    (
        maxPointField
        (
            error
            (
                cellToPoint(vFld),
                lowerRefineLevel,
                upperRefineLevel
            )
        )
    );

    // Mark cells that are candidates for refinement.
    forAll(cellError, cellI)
    {
        if (cellError[cellI] > 0)
        {
            candidateCell.set(cellI, 1);
        }
    }
}


labelList dynamicRefinePolyFvMesh::selectRefineCells
(
    const label maxCells,
    const label maxRefinement,
    const PackedBoolList& candidateCell
) const
{
    // Note: Assuming predominantly hex mesh, i.e. every refined cell causes 7
    // extra cells. This might lead to slight over shoot of maximum number of
    // cells.
    label nTotToRefine = (maxCells - globalData().nTotalCells()) / 7;

    const labelList& cellLevel = meshCutter_.cellLevel();

    // Count current selection
    label nCandidates = returnReduce(count(candidateCell, 1), sumOp<label>());

    // Collect all cells
    dynamicLabelList candidates(nCells());

    if (nCandidates < nTotToRefine)
    {
        // We won't exceed the maximum number of cells. Collect all candidate
        // cells that have refinement level smaller than max level
        forAll(candidateCell, cellI)
        {
            if
            (
                cellLevel[cellI] < maxRefinement
             && candidateCell.get(cellI) == 1
            )
            {
                candidates.append(cellI);
            }
        }
    }
    else
    {
        // We will exceed the maximum number of cells. Simply truncate the
        // list. It is possible to prefer certain cells based on error field.
        for (label level = 0; level < maxRefinement; level++)
        {
            forAll(candidateCell, cellI)
            {
                if
                (
                    cellLevel[cellI] == level
                 && candidateCell.get(cellI) == 1
                )
                {
                    candidates.append(cellI);
                }
            }

            if (returnReduce(candidates.size(), sumOp<label>()) > nTotToRefine)
            {
                break;
            }
        }
    }

    // Guarantee 2:1 refinement after refinement.

    // Create a label list for consistent set
    labelList consistentSet;

    // Note: the return type of consistentRefinement is Xfer<labelList>, so
    // we will transfer its contents into this list
    labelList transferredConsistentSet
    (
        meshCutter_.consistentRefinement
        (
            candidates.shrink(),
            true,                 // Add to set to guarantee 2:1,
            pointBasedRefinement_ // Whether to use point based refinement
        )
    );
    consistentSet.transfer(transferredConsistentSet);

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " cells for refinement out of " << globalData().nTotalCells()
        << "." << endl;

    return consistentSet;
}


labelList dynamicRefinePolyFvMesh::selectUnrefinePoints
(
    const scalar unrefineLevel,
    const PackedBoolList& markedCell,
    const scalarField& pFld
) const
{
    // All points that can be unrefined
    const labelList splitPoints(meshCutter_.getSplitPoints());

    dynamicLabelList newSplitPoints(splitPoints.size());

    forAll(splitPoints, i)
    {
        label pointI = splitPoints[i];

        if (pFld[pointI] < unrefineLevel)
        {
            // Check that all cells are not marked
            const labelList& pCells = pointCells()[pointI];

            bool hasMarked = false;

            forAll(pCells, pCellI)
            {
                if (markedCell.get(pCells[pCellI]) == 1)
                {
                    hasMarked = true;
                    break;
                }
            }

            if (!hasMarked)
            {
                newSplitPoints.append(pointI);
            }
        }
    }


    newSplitPoints.shrink();

    // Guarantee 2:1 refinement after unrefinement
    labelList consistentSet
    (
        meshCutter_.consistentUnrefinement
        (
            newSplitPoints,
            false,                // Remove from the set
            pointBasedRefinement_ // Whether to use point based unrefinement
        )
    );

    Info<< "Selected " << returnReduce(consistentSet.size(), sumOp<label>())
        << " split points out of a possible "
        << returnReduce(splitPoints.size(), sumOp<label>())
        << "." << endl;

    return consistentSet;
}


void dynamicRefinePolyFvMesh::extendMarkedCells(PackedBoolList& markedCell) const
{
    // Mark faces using any marked cell
    boolList markedFace(nFaces(), false);

    forAll(markedCell, cellI)
    {
        if (markedCell.get(cellI) == 1)
        {
            const cell& cFaces = cells()[cellI];

            forAll(cFaces, i)
            {
                markedFace[cFaces[i]] = true;
            }
        }
    }

    syncTools::syncFaceList(*this, markedFace, orEqOp<bool>(), false);

    // Update cells using any markedFace
    for (label faceI = 0; faceI < nInternalFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
            markedCell.set(faceNeighbour()[faceI], 1);
        }
    }
    for (label faceI = nInternalFaces(); faceI < nFaces(); faceI++)
    {
        if (markedFace[faceI])
        {
            markedCell.set(faceOwner()[faceI], 1);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

dynamicRefinePolyFvMesh::dynamicRefinePolyFvMesh(const IOobject& io)
:
    dynamicFvMesh(io),
    singleMotionUpdate_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
        .lookupOrDefault<Switch>("singleMotionUpdate", true)
    ),
    curTimeIndex_(-1),
    pointBasedRefinement_(false), // Set in the update() member function
    meshCutter_(*this),
    dumpLevel_(false),
    nRefinementIterations_(0)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

dynamicRefinePolyFvMesh::~dynamicRefinePolyFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool dynamicRefinePolyFvMesh::update()
{
    // Handling multiple calls in a single time step
    if
    (
        singleMotionUpdate_
     && curTimeIndex_ == this->time().timeIndex()
    )
    {
        // This is not the first call to update, simply return false
        return false;
    }

    // Update local time index
    curTimeIndex_ = this->time().timeIndex();

    // Re-read dictionary. Choosen since usually -small so trivial amount
    // of time compared to actual refinement. Also very useful to be able
    // to modify on-the-fly.
    dictionary refineDict
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE,
                false
            )
        ).subDict(typeName + "Coeffs")
    );

    label refineInterval = readLabel(refineDict.lookup("refineInterval"));

    bool hasChanged = false;

    if (refineInterval == 0)
    {
        changing(hasChanged);

        return false;
    }
    else if (refineInterval < 0)
    {
        FatalErrorIn("dynamicRefinePolyFvMesh::update()")
            << "Illegal refineInterval " << refineInterval << nl
            << "The refineInterval setting in the dynamicMeshDict should"
            << " be >= 1." << nl
            << exit(FatalError);
    }




    // Note: cannot refine at time 0 since no V0 present since mesh not
    //       moved yet.

    if (time().timeIndex() > 0 && time().timeIndex() % refineInterval == 0)
    {
        label maxCells = readLabel(refineDict.lookup("maxCells"));

        if (maxCells <= 0)
        {
            FatalErrorIn("dynamicRefinePolyFvMesh::update()")
                << "Illegal maximum number of cells " << maxCells << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        label maxRefinement = readLabel(refineDict.lookup("maxRefinement"));

        if (maxRefinement <= 0)
        {
            FatalErrorIn("dynamicRefinePolyFvMesh::update()")
                << "Illegal maximum refinement level " << maxRefinement << nl
                << "The maxCells setting in the dynamicMeshDict should"
                << " be > 0." << nl
                << exit(FatalError);
        }

        // Determining preferred consistency refinement based on maxRefinement
        if (maxRefinement > 2)
        {
            // Refinement level higher than 2 might produce 8:1 point
            // inconsistency, make sure that the point based consistency check
            // is used
            pointBasedRefinement_ = true;
        }

        // Check whether the user insists on specifying refinement strategy
        if (refineDict.found("pointBasedRefinement"))
        {
            const Switch userPointBasedRefinement =
                refineDict.lookup("pointBasedRefinement");

            // Check whether the user insisted on not using point based checking
            // and in case we determined that it should be used, issue a
            // warning and let the user have a try at it
            if (!userPointBasedRefinement && pointBasedRefinement_)
            {
                WarningIn("dynamicRefinePolyFvMesh::update()")
                    << "You are insisting on using face based consistency"
                    << " check for dynamic refinement."
                    << nl
                    << "Since you are allowing more than two maximum"
                    << " refinement levels, this might produce erroneous mesh"
                    << " due to 8:1 point conflicts."
                    << nl
                    << "In order to supress this message and use point based"
                    << " consistency checks, set pointBasedRefinement to on."
                    << endl;
            }

            pointBasedRefinement_ = userPointBasedRefinement;
        }

        const word fieldName(refineDict.lookup("field"));

        const volScalarField& vFld = lookupObject<volScalarField>(fieldName);

        const scalar lowerRefineLevel =
            readScalar(refineDict.lookup("lowerRefineLevel"));
        const scalar upperRefineLevel =
            readScalar(refineDict.lookup("upperRefineLevel"));
        const scalar unrefineLevel =
            readScalar(refineDict.lookup("unrefineLevel"));
        const label nBufferLayers =
            readLabel(refineDict.lookup("nBufferLayers"));

        // Cells marked for refinement
        PackedBoolList refineCell(nCells());

        if (globalData().nTotalCells() < maxCells)
        {
            // Determine candidates for refinement (looking at field only)
            selectRefineCandidates
            (
                lowerRefineLevel,
                upperRefineLevel,
                vFld,
                refineCell
            );

            // Select subset of candidates. Take into account max allowable
            // cells, refinement level, protected cells.
            labelList cellsToRefine
            (
                selectRefineCells
                (
                    maxCells,
                    maxRefinement,
                    refineCell
                )
            );

            label nCellsToRefine = returnReduce
            (
                cellsToRefine.size(), sumOp<label>()
            );

            if (nCellsToRefine > 0)
            {
                // Refine/update mesh and map fields
                autoPtr<mapPolyMesh> map = refine(cellsToRefine);

                // Update refineCell. Note that some of the marked ones have
                // not been refined due to constraints.
                {
                    const labelList& cellMap = map().cellMap();
                    const labelList& reverseCellMap = map().reverseCellMap();

                    PackedBoolList newRefineCell(cellMap.size());

                    forAll(cellMap, cellI)
                    {
                        label oldCellI = cellMap[cellI];

                        if (oldCellI < 0)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else if (reverseCellMap[oldCellI] != cellI)
                        {
                            newRefineCell.set(cellI, 1);
                        }
                        else
                        {
                            newRefineCell.set(cellI, refineCell.get(oldCellI));
                        }
                    }
                    refineCell.transfer(newRefineCell);
                }

                // Extend with a buffer layer to prevent neighbouring points
                // being unrefined.
                for (label i = 0; i < nBufferLayers; i++)
                {
                    extendMarkedCells(refineCell);
                }

                hasChanged = true;
            }
        }


        {
            // Select unrefineable points that are not marked in refineCell
            labelList pointsToUnrefine
            (
                selectUnrefinePoints
                (
                    unrefineLevel,
                    refineCell,
                    minCellField(vFld)
                )
            );

            label nSplitPoints = returnReduce
            (
                pointsToUnrefine.size(),
                sumOp<label>()
            );

            if (nSplitPoints > 0)
            {
                // Refine/update mesh
                unrefine(pointsToUnrefine);

                hasChanged = true;
            }
        }


        if ((nRefinementIterations_ % 10) == 0)
        {
            // Compact refinement history occassionally (how often?).
            // Unrefinement causes holes in the polyRefinementHistory.
            const_cast<polyRefinementHistory&>
            (
                meshCutter().history()
            ).compact();
        }
        nRefinementIterations_++;
    }

    changing(hasChanged);

    return hasChanged;
}


bool dynamicRefinePolyFvMesh::writeObject
(
    IOstream::streamFormat fmt,
    IOstream::versionNumber ver,
    IOstream::compressionType cmp
) const
{
    // Force refinement data to go to the current time directory.
    const_cast<polyRef&>(meshCutter_).setInstance(time().timeName());

    bool writeOk =
        dynamicFvMesh::writeObjects(fmt, ver, cmp)
     && meshCutter_.write();

    if (dumpLevel_)
    {
        volScalarField scalarCellLevel
        (
            IOobject
            (
                "cellLevel",
                time().timeName(),
                *this,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE,
                false
            ),
            *this,
            dimensionedScalar("level", dimless, 0)
        );

        const labelList& cellLevel = meshCutter_.cellLevel();

        forAll(cellLevel, cellI)
        {
            scalarCellLevel[cellI] = cellLevel[cellI];
        }

        writeOk = writeOk && scalarCellLevel.write();
    }

    return writeOk;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
