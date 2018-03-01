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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

\*---------------------------------------------------------------------------*/

#include "overlapFringe.H"
#include "oversetRegion.H"
#include "oversetMesh.H"
#include "polyPatchID.H"
#include "processorFvPatchFields.H"
#include "oversetFvPatchFields.H"
#include "typeInfo.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(overlapFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, overlapFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::overlapFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn("void overlapFringe::calcAddressing() const")
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get initial guess for holes and acceptors.
    // Algorithm:
    //    - Create indicator field for correct data exchange accross processor
    //      boundaries
    //    - Get holes from overset region and mark immediate neighbours of
    //      holes as acceptors
    //    - Loop through (optionally) user specified patches for
    //      initialising the overlap fringe assembly, marking face cells

    // Get necessary mesh data
    const fvMesh& mesh = region().mesh();
    const labelListList& cc = mesh.cellCells();

    // Note: Because we cannot assume anything about parallel decomposition and
    // we use neighbourhood walk algorithm, there is no easy way to go across
    // processor boundaries. We will create an indicator field marking all
    // possible cut cells and all possible face cells of given patches. Then, we
    // will use this indicator field to transfer the search to the other side.

    // Create the indicator field
    volScalarField processorIndicator
    (
        IOobject
        (
            "processorIndicator",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("minusOne", dimless, -1.0)
    );
    scalarField& processorIndicatorIn = processorIndicator.internalField();

    // Get cut holes from overset region
    const labelList& cutHoles = region().cutHoles();

    // Debug
    if (oversetMesh::debug && cutHoles.empty())
    {
        Pout<< "Did not find any holes to initialise the overlap fringe "
            << "assembly. Proceeding to patches..."
            << endl;
    }

    // Initialise mask field for eligible acceptors (cells that are not
    // holes)
    boolList eligibleAcceptors(mesh.nCells(), true);

    forAll (cutHoles, hI)
    {
        const label& holeCellI = cutHoles[hI];

        // Mask eligible acceptors
        eligibleAcceptors[holeCellI] = false;

        // Mark cut hole cell in processor indicator field
        processorIndicatorIn[holeCellI] = 1.0;
    }

    // Dynamic list for storing acceptors.
    // Note 1: capacity set to number of cells (trading off memory for
    // efficiency)
    // Note 2: inserting duplicates is avoided by updating eligibleAcceptors
    // mask
    dynamicLabelList candidateAcceptors(mesh.nCells());

    // Loop through cut holes and find acceptor candidates
    forAll (cutHoles, hI)
    {
        // Get neighbours of this hole cell
        const labelList& hNbrs = cc[cutHoles[hI]];

        // Loop through neighbours of this hole cell
        forAll (hNbrs, nbrI)
        {
            // Check whether the neighbouring cell is eligible
            const label& nbrCellI = hNbrs[nbrI];

            if (eligibleAcceptors[nbrCellI])
            {
                // Append the cell and mask it to avoid duplicate entries
                candidateAcceptors.append(nbrCellI);
                eligibleAcceptors[nbrCellI] = false;
            }
        }
    }

    // Debug
    if (oversetMesh::debug() && initPatchNames_.empty())
    {
        Pout<< "Did not find any specified patches to initialise the "
            << "overlap fringe assembly."
            << endl;
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Loop through patches and mark face cells as eligible acceptors
    forAll (initPatchNames_, nameI)
    {
        const polyPatchID curPatch
        (
            initPatchNames_[nameI],
            mesh.boundaryMesh()
        );

        if (!curPatch.active())
        {
            FatalErrorIn
            (
                "void overlapFringe::calcAddressing() const"
            )   << "Patch specified for fringe initialisation "
                << initPatchNames_[nameI] << " cannot be found"
                << abort(FatalError);
        }

        const unallocLabelList& curFaceCells =
            mesh.boundaryMesh()[curPatch.index()].faceCells();

        // Loop through face cells and mark candidate acceptors if
        // eligible
        forAll (curFaceCells, fcI)
        {
            // Get cell index
            const label& cellI = curFaceCells[fcI];

            // Mark acceptor face cell in processor indicator field
            processorIndicatorIn[cellI] = 1.0;

            // Check if the cell is eligible and if it is in region zone
            // (Note: the second check is costly)
            if
            (
                eligibleAcceptors[cellI]
             && rcz.whichCell(cellI) > -1
            )
            {
                candidateAcceptors.append(cellI);
                eligibleAcceptors[cellI] = false;
            }
        }
    }

    // Get boundary field
    volScalarField::GeometricBoundaryField& processorIndicatorBf =
        processorIndicator.boundaryField();

    // Perform update accross coupled boundaries, excluding overset patch
    evaluateNonOversetBoundaries(processorIndicatorBf);

    // Loop through boundary field
    forAll (processorIndicatorBf, patchI)
    {
        // Get patch field
        const fvPatchScalarField& chipf = processorIndicatorBf[patchI];

        // Only perform acceptor search if this is a processor boundary
        if (isA<processorFvPatchScalarField>(chipf))
        {
            // Get neighbour field
            const scalarField nbrProcIndicator =
                chipf.patchNeighbourField();

            // Get face cells
            const unallocLabelList& fc = chipf.patch().faceCells();

            // Loop through neighbouring processor field
            forAll (nbrProcIndicator, pfaceI)
            {
                if
                (
                    nbrProcIndicator[pfaceI] > 0.0
                 && eligibleAcceptors[fc[pfaceI]]
                )
                {
                    // The cell on the other side is a hole or acceptor from
                    // face cells of a given patch, while the cell on this side
                    // has not been marked yet neither as an acceptor or as a
                    // hole. Append the cell to candidate acceptors and mark it
                    // as ineligible in order to propage the fringe on this side
                    candidateAcceptors.append(fc[pfaceI]);
                    eligibleAcceptors[fc[pfaceI]] = false;
                }
            }
        }
    }

    // Now we have a decent first guess for acceptors that will be used as
    // an initial condition for the iterative overlap assembly
    // process.
    // Transfer the acceptor list and allocate empty fringeHoles list, which
    // may be populated in updateIteration member function
    acceptorsPtr_ = new labelList(candidateAcceptors.xfer());
    fringeHolesPtr_ = new labelList(cutHoles);
}


void Foam::overlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
}


void Foam::overlapFringe::evaluateNonOversetBoundaries
(
    volScalarField::GeometricBoundaryField& psib
) const
{
    // Code practically copy/pasted from GeometricBoundaryField::evaluateCoupled
    // GeometricBoundaryField should be redesigned to accomodate for such needs
    if
    (
        Pstream::defaultComms() == Pstream::blocking
     || Pstream::defaultComms() == Pstream::nonBlocking
    )
    {
        forAll(psib, patchI)
        {
            // Get fvPatchField
            fvPatchScalarField& psip = psib[patchI];

            if (psip.coupled() && !isA<oversetFvPatchScalarField>(psip))
            {
                psip.initEvaluate(Pstream::defaultComms());
            }
        }

        // Block for any outstanding requests
        if (Pstream::defaultComms() == Pstream::nonBlocking)
        {
            IPstream::waitRequests();
            OPstream::waitRequests();
        }

        forAll(psib, patchI)
        {
            // Get fvPatchField
            fvPatchScalarField& psip = psib[patchI];

            if (psip.coupled() && !isA<oversetFvPatchScalarField>(psip))
            {
                psip.evaluate(Pstream::defaultComms());
            }
        }
    }
    else if (Pstream::defaultComms() == Pstream::scheduled)
    {
        const lduSchedule& patchSchedule =
            region().mesh().globalData().patchSchedule();

        forAll(patchSchedule, patchEvalI)
        {
            if (patchSchedule[patchEvalI].init)
            {
                // Get fvPatchField
                fvPatchScalarField psip = psib[patchSchedule[patchEvalI].patch];

                if (psip.coupled() && !isA<oversetFvPatchScalarField>(psip))
                {
                    psip.initEvaluate(Pstream::scheduled);
                }
            }
            else
            {
                // Get fvPatchField
                fvPatchScalarField psip = psib[patchSchedule[patchEvalI].patch];

                if (psip.coupled() && !isA<oversetFvPatchScalarField>(psip))
                {
                    psip.evaluate(Pstream::scheduled);
                }
            }
        }
    }
    else
    {
        FatalErrorIn("overlapFringe::evaluateNonOversetBoundaries()")
            << "Unsuported communications type "
            << Pstream::commsTypeNames[Pstream::defaultCommsType()]
            << exit(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::overlapFringe::overlapFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    fringeHolesPtr_(NULL),
    acceptorsPtr_(NULL),
    finalDonorAcceptorsPtr_(NULL),

    donorSuitability_
    (
        donorSuitability::donorSuitability::New(*this, dict)
    ),
    initPatchNames_
    (
        dict.lookupOrDefault<wordList>("initPatchNames", wordList())
    ),
    minGlobalFraction_
    (
        readScalar(dict.lookup("suitablePairFraction"))
    ),
    cumulativeDonorAcceptorsPtr_(NULL),
    cacheFringe_(dict.lookupOrDefault<Switch>("cacheFringe", false)),
    fringeIter_(0)
{
    // Sanity check
    if (minGlobalFraction_ < SMALL || minGlobalFraction_ > 1.0)
    {
        FatalErrorIn
        (
            "overlapFringe::overlapFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict,\n"
            ")\n"
        )   << "Invalid suitablePairFraction found while reading the overlap "
            << "fringe dictionary."
            << nl
            << "Please specify value between 0 and 1."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::overlapFringe::~overlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::overlapFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    if (!fringeHolesPtr_ || !acceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::updateIteration()")
            << "fringeHolesPtr_ or acceptorsPtr_ is not allocated. "
            << "Make sure you have called acceptors() or fringeHoles() to "
            << "calculate the initial set of donor/acceptors before "
            << "actually updating iteration."
            << abort(FatalError);
    }

    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::updateIteration()")
            << "Called iteration update with finalDonorAcceptorsPtr_ "
            << "allocated. This means that the final overlap has been "
            << "achieved, prohibiting calls to updateIteration."
            << abort(FatalError);
    }

    // Increment iteration counter for output
    ++fringeIter_;

    // Allocate worker cumulative donor/acceptor list if it has not been
    // allocated yet (first iteration). Use largest possible size to prevent
    // any resizing
    if (!cumulativeDonorAcceptorsPtr_)
    {
        cumulativeDonorAcceptorsPtr_ = new donorAcceptorDynamicList
        (
            region().mesh().nCells()
        );
    }
    donorAcceptorDynamicList& cumDAPairs = *cumulativeDonorAcceptorsPtr_;

    // Create a list containing unsuitable donors
    donorAcceptorDynamicList unsuitableDAPairs(donorAcceptorRegionData.size());

    // Loop through donor/acceptor pairs and perform mark-up
    forAll (donorAcceptorRegionData, daPairI)
    {
        if
        (
            donorSuitability_->isDonorSuitable(donorAcceptorRegionData[daPairI])
        )
        {
            // Donor is suitable, add it directly to the cumulative list
            cumDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
        else
        {
            // Donor is not suitable, append it to the unsuitable list
            unsuitableDAPairs.append(donorAcceptorRegionData[daPairI]);
        }
    }

    // Calculate the number of total suitable pairs found so far and the number
    // of total pairs
    const label nSuitablePairs =
        returnReduce<label>(cumDAPairs.size(), sumOp<label>());

    const label nTotalPairs = nSuitablePairs
      + returnReduce<label>(unsuitableDAPairs.size(), sumOp<label>());

    const scalar suitabilityFrac = scalar(nSuitablePairs)/scalar(nTotalPairs);

    // Print information
    Info<< "Overlap fringe iteration: " << fringeIter_
        << " for region: " << region().name()
        << nl
        << "Cumulative suitable pairs: " << nSuitablePairs
        << ", total number of pairs: " << nTotalPairs
        << " (" << suitabilityFrac*100 << "%)"
        << endl;

    // Check whether the criterion has been satisfied
    if (suitabilityFrac > minGlobalFraction_)
    {
        // Append unsuitable donors to the list as well
        cumDAPairs.append(unsuitableDAPairs);

        // Now that we have reached suitability criterion specified by the user,
        // we need to clean up a bit. Namely, it is possible that a certain
        // acceptor cell is completely surrounded by holes or other acceptor, so
        // this cell needs to become a hole as well. For easier parallel
        // processing, we will create an indicator field where hole and acceptor
        // cells are marked with 1 and all the other cells (live cells) are
        // marked with -1. We will then use this indicator field to determine
        // whether this acceptor needs to become a hole.

        // Get mesh
        const fvMesh& mesh = region().mesh();

        // Create the processor indicator field to transfer hole cells to the
        // other side
        volScalarField holeIndicator
        (
            IOobject
            (
                "holeIndicator_" + region().name(),
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("minusOne", dimless, -1.0)
        );
        scalarField& holeIndicatorIn = holeIndicator.internalField();

        // Transfer fringeHolesPtr into the dynamic list for efficiency. Note:
        // will be transfered back at the end of the scope.
        dynamicLabelList allFringeHoles(fringeHolesPtr_->xfer());

        // Loop through all fringe holes and mark them
        forAll (allFringeHoles, hcI)
        {
            holeIndicatorIn[allFringeHoles[hcI]] = 1.0;
        }

        // Loop through all acceptors and mark them
        forAll (cumDAPairs, daPairI)
        {
            holeIndicatorIn[cumDAPairs[daPairI].acceptorCell()] = 1.0;
        }

        // Get boundary field
        volScalarField::GeometricBoundaryField& holeIndicatorb =
            holeIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        evaluateNonOversetBoundaries(holeIndicatorb);

        // Get necessary mesh data
        const cellList& meshCells = mesh.cells();
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // List of acceptors to be converted to holes
        boolList accBecomingHoles(cumDAPairs.size(), false);

        // Loop through all donor/acceptor pairs collected so far
        forAll (cumDAPairs, daPairI)
        {
            // Get acceptor cell index
            const label& accI = cumDAPairs[daPairI].acceptorCell();

            // Get faces of this cell
            const cell& accFaces = meshCells[accI];

            // Create a bool whether this acceptor needs to be converted to hole
            bool convertToHole = true;

            // Loop through faces
            forAll (accFaces, faceI)
            {
                // Get global face index
                const label& gfI = accFaces[faceI];

                // Check whether this is an internal face or patch face
                if (mesh.isInternalFace(gfI))
                {
                    // Internal face, check whether I'm owner or neighbour
                    if (own[gfI] == accI)
                    {
                        // I'm owner, check whether the neighbour is live
                        if (holeIndicatorIn[nei[gfI]] < 0.0)
                        {
                            // This acceptor has a live cell for neighbour,
                            // update the flag and continue
                            convertToHole = false;
                            continue;
                        }
                    }
                    else
                    {
                        // I'm neighbour, check whether the owner is live
                        if (holeIndicatorIn[own[gfI]] < 0.0)
                        {
                            // This acceptor has a live cell for neighbour,
                            // update the flag and continue
                            convertToHole = false;
                            continue;
                        }
                    }
                }
                else
                {
                    // Get patch and face index
                    const label patchI = mesh.boundaryMesh().whichPatch(gfI);
                    const label pfI =
                        mesh.boundaryMesh()[patchI].whichFace(gfI);

                    // Only consider processor patches
                    if (isA<processorPolyPatch>(mesh.boundaryMesh()[patchI]))
                    {
                        // Note: patch stores neighbour field after evaluation
                        if (holeIndicatorb[patchI][pfI] < 0.0)
                        {
                            // This acceptor has a live cell for neighbour on
                            // the other processor, update the flag and continue
                            convertToHole = false;
                            continue;
                        }
                    }
                }
            }

            // Mark whether this acceptor cell has to be converted to hole
            accBecomingHoles[daPairI] = convertToHole;
        }

        // Now we need to filter the data: append acceptors that need to be
        // converted to holes into allFringeHoles and insert all other acceptors
        // into finalDAPairs temporary container
        // Create another dynamic list to collect final donor/acceptor pairs
        donorAcceptorDynamicList finalDAPairs(cumDAPairs.size());

        // Count number of acceptor holes that need to be converted to holes
        label nAccToHoles = 0;

        // Loop all current donor/acceptor pairs
        forAll(cumDAPairs, daPairI)
        {
            if (accBecomingHoles[daPairI])
            {
                // Append the acceptor to list of holes
                allFringeHoles.append(cumDAPairs[daPairI].acceptorCell());
                ++nAccToHoles;
            }
            else
            {
                // Append the donor/acceptor pair to finalDAPairs list
                finalDAPairs.append(cumDAPairs[daPairI]);
            }
        }

        // Tranfer back the allFringeHoles dynamic list into member data
        fringeHolesPtr_->transfer(allFringeHoles);

        // Transfer ownership of the final donor/acceptor list to the
        // finalDonorAcceptorsPtr_
        finalDonorAcceptorsPtr_ = new donorAcceptorList
        (
            finalDAPairs.xfer()
        );

        // At least 100*minGlobalFraction_ % of suitable donor/acceptor pairs
        // have been found.
        Info<< "Converted " << nAccToHoles << " acceptors to holes."
            << nl
            << "Finished assembling overlap fringe. " << endl;

        // Set the flag to true
        updateSuitableOverlapFlag(true);
    }
    else
    {
        // A sufficient number of suitable donor/acceptors has not been
        // found. Go through unsuitable donor/acceptor pairs and find a new
        // batch of acceptors and holes for the next iteration

        // Get necessary mesh data
        const fvMesh& mesh = region().mesh();
        const labelListList& cc = mesh.cellCells();

        // Create the processor indicator field to transfer the unsuitable
        // acceptors to the other side
        volScalarField processorIndicator
        (
            IOobject
            (
                "processorIndicator",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar("minusOne", dimless, -1.0)
        );
        scalarField& processorIndicatorIn = processorIndicator.internalField();

        // Transfer fringeHolesPtr into the dynamic list for efficiency. Note:
        // will be transfered back at the end of the scope.
        dynamicLabelList cumFringeHoles(fringeHolesPtr_->xfer());

        // Create mask to prevent wrong and duplicate entries (i.e. we cannot
        // search backwards through existing acceptors and holes)
        boolList freeCells(mesh.nCells(), true);

        // Mask all considered suitable acceptor cells so far
        forAll (cumDAPairs, cpI)
        {
            freeCells[cumDAPairs[cpI].acceptorCell()] = false;
        }

        // Mask all current unsuitable acceptor pairs as well
        forAll (unsuitableDAPairs, upI)
        {
            const label& accCellI = unsuitableDAPairs[upI].acceptorCell();

            freeCells[accCellI] = false;

            // Mark unsuitable pair for possible processor transfer
            processorIndicatorIn[accCellI] = 1.0;
        }

        // Mask all fringe holes
        forAll (cumFringeHoles, cfhI)
        {
            const label& fhCellI = cumFringeHoles[cfhI];

            freeCells[fhCellI] = false;

            // Mark fringe hole for possible processor transfer
            processorIndicatorIn[fhCellI] = 1.0;
        }

        // Create dynamic list to efficiently append new batch of
        // acceptors. Note: allocate enough storage.
        dynamicLabelList newAcceptors(10*unsuitableDAPairs.size());

        // Loop through unsuitable acceptors
        forAll (unsuitableDAPairs, upI)
        {
            // Get acceptor cell and its neighbours
            const label& accI = unsuitableDAPairs[upI].acceptorCell();
            const labelList& aNbrs = cc[accI];

            // Loop through neighbours of this acceptor cell
            forAll (aNbrs, nbrI)
            {
                // Check whether the neighbouring cell is free
                const label& nbrCellI = aNbrs[nbrI];

                if (freeCells[nbrCellI])
                {
                    // This cell is neither an old acceptor, fringe hole nor it
                    // has been considered previously. Append it to the
                    // newAcceptors list and mark it as visited
                    newAcceptors.append(nbrCellI);
                    freeCells[nbrCellI] = false;
                }
            }

            // Append this "old" acceptor cell into fringe holes list
            cumFringeHoles.append(accI);
        }

        // Transfer the fringe accross possible processor boundaries

        // Get boundary field
        volScalarField::GeometricBoundaryField& processorIndicatorBf =
            processorIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        evaluateNonOversetBoundaries(processorIndicatorBf);

        // Loop through boundary field
        forAll (processorIndicatorBf, patchI)
        {
            // Get patch field
            const fvPatchScalarField& chipf = processorIndicatorBf[patchI];

            // Only perform acceptor search if this is a processor boundary
            if (isA<processorFvPatchScalarField>(chipf))
            {
                // Get neighbour field
                const scalarField nbrProcIndicator =
                    chipf.patchNeighbourField();

                // Get face cells
                const unallocLabelList& fc = chipf.patch().faceCells();

                // Loop through neighbouring processor field
                forAll (nbrProcIndicator, pfaceI)
                {
                    if
                    (
                        nbrProcIndicator[pfaceI] > 0.0
                     && freeCells[fc[pfaceI]]
                    )
                    {
                        // The cell on the other side is a hole or acceptor,
                        // while the cell on this side has not been marked yet
                        // as an acceptor. Append the cell to new set of
                        // acceptors and mark it as ineligible in order to
                        // propage the fringe on this side
                        newAcceptors.append(fc[pfaceI]);
                        freeCells[fc[pfaceI]] = false;
                    }
                }
            }
        }

        // Transfer back cumulative fringe holes into the fringeHolesPtr_
        fringeHolesPtr_->transfer(cumFringeHoles);

        // Transfer new acceptors into the acceptors list
        acceptorsPtr_->transfer(newAcceptors);

        // Set the flag to false (suitable overlap not found)
        updateSuitableOverlapFlag(false);
    }

    return foundSuitableOverlap();
}



const Foam::labelList& Foam::overlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    // Debug, write fringe holes as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor fringe holes into a cell set." << endl;

        cellSet holesSet
        (
            mesh(),
            "fringeHolesProc" + Pstream::myProcNo() + region().name(),
            labelHashSet(*fringeHolesPtr_)
        );

        holesSet.write();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::overlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    // Debug, write candidate acceptors as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor candidate acceptors into a cell set." << endl;

        cellSet candidateAcceptorsSet
        (
            mesh(),
            "candidateAcceptorsProc" + Pstream::myProcNo() + region().name(),
            labelHashSet(*acceptorsPtr_)
        );

        candidateAcceptorsSet.write();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::overlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("overlapFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called overlapFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("overlapFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    if (cacheFringe_)
    {
        // Get reference to final donor/acceptor pairs
        const donorAcceptorList& finalDAPairs = *finalDonorAcceptorsPtr_;

        // Clear acceptors
        deleteDemandDrivenData(acceptorsPtr_);

        // Now allocate with the correct size. Note: since it is expected that
        // acceptorsPtr_->size() (before destruction) is smaller than
        // finalDonorAcceptorsPtr_.size(), this destruction and initialization
        // should not represent an overhead.
        acceptorsPtr_ = new labelList(finalDAPairs.size());
        labelList& acceptors = *acceptorsPtr_;

        // Set acceptors for the next fringe assembly process
        forAll (finalDAPairs, daPairI)
        {
            acceptors[daPairI] = finalDAPairs[daPairI].acceptorCell();
        }

        // Note: fringe holes actually hold the complete list, there's nothing
        // to do
        Info<< "Cached "
            << returnReduce<label>(acceptorsPtr_->size(), sumOp<label>())
            << " acceptors and "
            << returnReduce<label>(fringeHolesPtr_->size(), sumOp<label>())
            << " fringe holes for region: " << region().name()
            << endl;
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::overlapFringe::update() const
{
    Info<< "overlapFringe::update() const" << endl;

    if (cacheFringe_)
    {
        // If the cache is switched on, simply do not clear acceptorsPtr_ and
        // fringeHolesPtr_ which now hold all final acceptors and holes that
        // will be used to start the next iteration

        // Now clear final and cumulative donor acceptors
        deleteDemandDrivenData(finalDonorAcceptorsPtr_);
        deleteDemandDrivenData(cumulativeDonorAcceptorsPtr_);
    }
    else
    {
        // Clear everything, including acceptors and fringe holes
        clearAddressing();
    }

    // Reset iteration counter
    fringeIter_ = 0;

    // Set flag to false
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
