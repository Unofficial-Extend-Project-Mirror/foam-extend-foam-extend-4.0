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

#include "adaptiveOverlapFringe.H"
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
    defineTypeNameAndDebug(adaptiveOverlapFringe, 0);
    addToRunTimeSelectionTable
    (
        oversetFringe,
        adaptiveOverlapFringe,
        dictionary
    );
}


// * * * * * * * * * * * * * Static Member Functions  * * * * * * * * * * * * //

bool Foam::adaptiveOverlapFringe::sortBySuitabilityValue
(
    const iterationData& lhs,
    const iterationData& rhs
)
{
    return lhs.suitability() < rhs.suitability();
}


void Foam::adaptiveOverlapFringe::suitabilityFractionSlope
(
    FIFOStack<iterationData> iterHist,
    scalar& alpha,
    scalar& beta
) const
{
    // Linear regression coefficients y = alpha + beta*x
    // beta = sum((x_i - x_mean)*(y_i - y_mean))/sum(x_i - x_mean)^2
    // alpha = y_mean - beta*x_mean
    // x -> iteration number
    // y -> suitability

    // Calculate mean values first
    scalar iterMean = 0;
    scalar suitabilityMean = 0;

    forAllConstIter(FIFOStack<iterationData>, iterHist, it)
    {
        iterMean += it().iteration();
        suitabilityMean += it().suitability();
    }

    reduce(iterMean, sumOp<scalar>());
    iterMean /= iterHist.size();
    reduce(suitabilityMean, sumOp<scalar>());
    suitabilityMean /= iterHist.size();

    // Calculate denominator and numerator for beta
    scalar n = 0;
    scalar d = 0;

    forAllConstIter(FIFOStack<iterationData>, iterHist, it)
    {
        n += (it().iteration() - iterMean)*
             (it().suitability() - suitabilityMean);
        d += Foam::sqr((it().iteration() - iterMean));
    }

    reduce(n, sumOp<scalar>());
    reduce(d, sumOp<scalar>());

    // Calculate regression coefficients
    beta = n/d;
    alpha = suitabilityMean - beta*iterMean;
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::adaptiveOverlapFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn("void adaptiveOverlapFringe::calcAddressing() const")
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Get initial guess for holes and acceptors.
    // Algorithm:
    //    - Create indicator field for correct data exchange accross processor
    //      boundaries
    //    - Get holes from overset region (and optionally from specified set)
    //      and mark immediate neighbours of holes as acceptors
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

    // Read user specified holes into allHoles list. Note: use cellZone rather
    // than cellSet to have correct behaviour on dynamic mesh simulations
    // We will silently proceed if the zone is not found since this option is
    // not mandatory but is useful in certain cases

    // Get zone index
    const label zoneID = mesh.cellZones().findZoneID(holesZoneName_);

    // Create a hash set for allHoles
    labelHashSet allHoles;

    if (zoneID > -1)
    {
        // Get the zone for holes and append them to set
        const labelList& specifiedHoles = mesh.cellZones()[zoneID];

        allHoles.insert(specifiedHoles);
    }
    // else silently proceed without user-specified holes

    // Extend allHoles with cutHoles
    forAll (cutHoles, chI)
    {
        // Note: duplicated are removed because we're using hash set
        allHoles.insert(cutHoles[chI]);
    }

    // Mark all holes
    forAllConstIter (labelHashSet, allHoles, iter)
    {
        const label& holeCellI = iter.key();

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

    // Loop through all holes and find acceptor candidates
    forAllConstIter (labelHashSet, allHoles, iter)
    {
        // Get neighbours of this hole cell
        const labelList& hNbrs = cc[iter.key()];

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
                "void adaptiveOverlapFringe::calcAddressing() const"
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
    overlapFringe::evaluateNonOversetBoundaries(processorIndicatorBf);

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

    // Issue an error if no acceptors have been found for initial guess
    if (returnReduce(candidateAcceptors.size(), sumOp<label>()) == 0)
    {
        FatalErrorIn
        (
            "void adaptiveOverlapFringe::calcAddressing() const"
        )   << "Did not find any acceptors to begin with."
            << "Check definition of adaptiveOverlap in oversetMeshDict"
            << " for region: " << this->region().name() << nl
            << "More specifically, check definition of:" << nl
            << "1. holePatches (mandatory entry)" << nl
            << "2. holes (optional entry)" << nl
            << "3. initPatchNames (optional entry)"
            << abort(FatalError);
    }


    // Now we have a decent first guess for acceptors that will be used as
    // an initial condition for the iterative overlap assembly
    // process.
    // Transfer the acceptor list and allocate empty fringeHoles list, which
    // may be populated in updateIteration member function
    acceptorsPtr_ = new labelList(candidateAcceptors.xfer());
    fringeHolesPtr_ = new labelList(allHoles.toc().xfer());
}


void Foam::adaptiveOverlapFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);

    // Reset and clear out helper data
    suitableDAPairs_.clear();
    suitablePairsSuit_ = 0;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::adaptiveOverlapFringe::adaptiveOverlapFringe
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

    holesZoneName_(dict.lookupOrDefault<word>("holes", word())),
    initPatchNames_
    (
        dict.lookupOrDefault<wordList>("initPatchNames", wordList())
    ),

    donorSuitability_
    (
        donorSuitability::donorSuitability::New(*this, dict)
    ),
    fringeIter_(0),
    specifiedIterationsNumber_
    (
        dict.lookupOrDefault<label>("specifiedIterationsNumber", 5)
    ),
    minSuitabilityRate_
    (
        dict.lookupOrDefault<scalar>("minSuitabilityRate", 0.0)
    ),
    maxIter_
    (
        dict.lookupOrDefault<label>("maximumIterations", 100)
    ),
    relativeCounter_(specifiedIterationsNumber_),
    additionalIterations_
    (
        dict.lookupOrDefault<Switch>("additionalIterations", true)
    ),
    orphanSuitability_
    (
        dict.lookupOrDefault<scalar>("orphanSuitability", 1)
    ),
    suitablePairsSuit_(0)

{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::adaptiveOverlapFringe::~adaptiveOverlapFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::adaptiveOverlapFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    if (!fringeHolesPtr_ || !acceptorsPtr_)
    {
        FatalErrorIn("adaptiveOverlapFringe::updateIteration()")
            << "fringeHolesPtr_ or acceptorsPtr_ is not allocated. "
            << "Make sure you have called acceptors() or fringeHoles() to "
            << "calculate the initial set of donor/acceptors before "
            << "actually updating iteration."
            << abort(FatalError);
    }

    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("adaptiveOverlapFringe::updateIteration()")
            << "Called iteration update with finalDonorAcceptorsPtr_ "
            << "allocated. This means that the final overlap has been "
            << "achieved, prohibiting calls to updateIteration."
            << abort(FatalError);
    }

    // Increment iteration counter for output
    ++fringeIter_;

    // Print info
    Info<< "Region: " << region().name() << ", iteration: "
        << fringeIter_ << endl;

    // Store donor/acceptor pairs whose donors are not within bounding box or
    // their suitability is lower than threshold into unsuitableDAPairs
    // donorAcceptorDynamicList. Neighbours of those acceptors will be
    // candidates for new acceptors.
    donorAcceptorDynamicList unsuitableDAPairs(donorAcceptorRegionData.size());

    // All donor/acceptor pairs from one iteration. Store this list into
    // iterationDataHistory_ after average suitability is calculated.
    donorAcceptorDynamicList storageDAPairs(donorAcceptorRegionData.size());

    // Suitability fraction variables

    // Not within bounding box pairs counter
    label notWithinBBCounter = 0;

    // Unsuitable pairs cumulative suitability
    scalar unsuitableSuitCum = 0;

    // Loop through donor/acceptor pairs and divide received donor/acceptor
    // pairs into two lists:
    //     1) unsuitableDAPairs,
    //     2) suitableDAPairs_ - those pairs are valid, keep them
    //        through the whole iterative process.
    forAll (donorAcceptorRegionData, daPairI)
    {
        // Get current donor/acceptor pair
        const donorAcceptor& curDA = donorAcceptorRegionData[daPairI];

        if (!curDA.withinBB())
        {
            // Donor of this acceptor is not within bounding box.
            // Append this pair to unsuitableDAPairs list.
            unsuitableDAPairs.append(curDA);

            // Add suitability of this pair to cumulative suitability.
            unsuitableSuitCum -= orphanSuitability_;

            // Increment number of pairs
            ++notWithinBBCounter;
        }
        else
        {
            // Those donor/acceptor pairs are valid, i.e. donor is within
            // bounding box

            // Calculate donor acceptor suitability
            const scalar donorAcceptorSuit =
                donorSuitability_->suitabilityFraction(curDA);

            if (!donorSuitability_->isDonorSuitable(curDA))
            {
                // Suitability of this pair is lower than threshold.
                // Append it to unsuitableDAPairs list.
                unsuitableDAPairs.append(curDA);

                // For debug
                unsuitableSuitCum += donorAcceptorSuit;
            }
            else
            {
                // Suitability of this pair is greater than threshold.

                // Add suitability to suitablePairsSuit_
                suitablePairsSuit_ += donorAcceptorSuit;

                // Append pair to suitableDAPairs_
                suitableDAPairs_.append(curDA);
            }
        }
    }

    // Append suitable and unsuitable pairs to storageDAPairs
    storageDAPairs.append(unsuitableDAPairs);
    storageDAPairs.append(suitableDAPairs_);

    // Reduce all necessary information
    reduce(notWithinBBCounter, sumOp<label>());
    reduce(unsuitableSuitCum, sumOp<scalar>());

    const scalar globalSuitablePairsSuit =
        returnReduce(suitablePairsSuit_, sumOp<scalar>());
    const label nGlobalDAPairs =
        returnReduce(storageDAPairs.size(), sumOp<label>());
    const label nGlobalUnsuitablePairs =
        returnReduce(unsuitableDAPairs.size(), sumOp<label>());

    // Calculate donor/acceptor average suitability from current iteration
    const scalar donorAcceptorSuitAverage =
        (unsuitableSuitCum + globalSuitablePairsSuit)/nGlobalDAPairs;

    // Copy fringe holes list in order to store it into iterationDataObject
    labelList allFringeHoles(*fringeHolesPtr_);

    // Store everything that is needed to choose and reconstruct the best
    // overlap fringe assembly into iterationDataObject and push into stack:
    //      1) donor/acceptor pairs,
    //      2) fringe holes,
    //      3) average suitability,
    //      4) iteration ordinal number.
    iterationDataHistory_.push
    (
        iterationData
        (
            storageDAPairs.xfer(),
            allFringeHoles.xfer(),
            donorAcceptorSuitAverage,
            fringeIter_
        )
    );

    Info<< "Average donor/acceptor pair suitability: "
        << donorAcceptorSuitAverage*100 << " %."
        << endl;

    // Debug: write stored acceptors and holes and some general information
    if (adaptiveOverlapFringe::debug)
    {
        // Print information
        Info<< "Found " << notWithinBBCounter << " pairs that are not within "
            << "bounding box and "
            << nGlobalUnsuitablePairs - notWithinBBCounter
            << " pairs whose suitability " << nl
            << "is lower than minLocalSuit. "
            << endl;

        // Iteration ordinal number from max object
        label iterNum = fringeIter_;

        // Acceptors field from max object
        volScalarField storedAcc
        (
            IOobject
            (
                "AccStor" + Foam::name(iterNum) + region().name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1)
        );
        scalarField& acceptorsIn = storedAcc.internalField();

        // Holes field from max object
        volScalarField holes
        (
            IOobject
            (
                "Holes" + Foam::name(iterNum) + region().name(),
                mesh().time().timeName(),
                mesh(),
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh(),
            dimensionedScalar("one", dimless, 1)
        );
        scalarField& holesIn = holes.internalField();

        // Mark all valid acceptors
        forAll(suitableDAPairs_, aI)
        {
            acceptorsIn[suitableDAPairs_[aI].acceptorCell()] = 2;
        }

        // Mark all unsuitable acceptors
        forAll(unsuitableDAPairs, aI)
        {
            acceptorsIn[unsuitableDAPairs[aI].acceptorCell()] = 0;
        }

        // Mark all holes
        forAll(allFringeHoles, hI)
        {
            // Get a hole cell index
            const label holeIndex = allFringeHoles[hI];

            // Mark hole
            holesIn[holeIndex] = -1;
        }

        // Write fields
        storedAcc.write();
        holes.write();
    }

    // Calculate slope if:
    //      1) User specified that he wants additional iterations and
    //      2) the number of made iterations is equal to the number of specified
    //         iterations and
    //      3) specified iterations number is greater than 1.
    // If the number of specified iterations is 1, I assume that user wants
    // only 1 iteration.
    // Note: To calculate linear regression coefficients, i.e. slope and
    // intercept, we need at least two points. So, calculating linear
    // regression coefficients with 1 iteration would lead to floating
    // point exception.

    if
    (
        (additionalIterations_)
     && (fringeIter_ == relativeCounter_)
     && (specifiedIterationsNumber_ > 1)
     && (fringeIter_ < maxIter_)
    )
    {
        // Slope (i.e. if y = alpha + beta*x, slope is beta)
        scalar beta = 0;

        // Intercept (i.e. if y = alpha + beta*x, intercept is alpha)
        scalar alpha = 0;

        // Calculate slope and intercept
        suitabilityFractionSlope
        (
            iterationDataHistory_,
            alpha,
            beta
        );

        // Print information
        if (adaptiveOverlapFringe::debug)
        {
            Info<< "Average donor suitability fraction depending on"
                << " iteration number: " << nl
                << "ADSF = " << alpha*100 << " + " << beta*100 << "*iter."
                << nl << endl;
        }


        // If the donor/acceptor suitability gradient is positive,
        // make 1 additional iteration
        if (beta > minSuitabilityRate_)
        {
            // Pop the bottom element of the stack
            iterationDataHistory_.pop();

            // Increment relative counter
            ++relativeCounter_;

            // Print information
            if (adaptiveOverlapFringe::debug)
            {
                Info<< "1 additional iteration for region " << region().name()
                    << " will be performed." << endl;
            }
        }
    }

    // Go through unsuitable donor/acceptor pairs from previous iteration and
    // find a new batch of acceptors and holes for the next iteration if:
    // 1) Desired number of iterations has not been made or
    //    relativeCounter_ is incremented because suitability gradient was
    //    positive and
    // 2) There is a number of  unsuitable donor/acceptor pairs what
    //    indicates that suitable overlap is found.
    if
    (
        (fringeIter_ < relativeCounter_)
     && (nGlobalUnsuitablePairs != 0)
    )
    {
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

        // Mask all current unsuitable acceptor pairs
        forAll (unsuitableDAPairs, upI)
        {
            const label& accCellI = unsuitableDAPairs[upI].acceptorCell();

            freeCells[accCellI] = false;

            // Mark unsuitable pair for possible processor transfer
            processorIndicatorIn[accCellI] = 1.0;
        }

        // Mask all current suitable acceptor pairs
        forAll (suitableDAPairs_, upI)
        {
            const label& accCellI = suitableDAPairs_[upI].acceptorCell();

            freeCells[accCellI] = false;

            // Mark suitable pair for possible processor transfer
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
        dynamicLabelList newAcceptors(2*donorAcceptorRegionData.size());

        // Loop through unsuitable donor/acceptor pairs.
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
                    // This cell is neither an old acceptor, fringe hole
                    // nor it has been considered previously. Append it to
                    // the newAcceptors list and mark it as visited
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
        overlapFringe::evaluateNonOversetBoundaries(processorIndicatorBf);

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

        if (returnReduce(newAcceptors.empty(), andOp<bool>()))
        {
            FatalErrorIn("adaptiveOverlapFringe::updateIteration()")
                << "Did not find any new candidate acceptors."
                << nl
                << "Please review your overlap fringe assembly settings."
                << abort(FatalError);
        }

        // Transfer back cumulative fringe holes into the fringeHolesPtr_
        fringeHolesPtr_->transfer(cumFringeHoles);

        // Transfer new acceptors into the acceptors list
        acceptorsPtr_->transfer(newAcceptors);

        // Set the flag to false (suitable overlap not found)
        updateSuitableOverlapFlag(false);
    }
    else
    {
        // Specified number of iterations is made or suitable overlap is found.
        // There is no possibility for memory leak because
        // Foam::adaptiveOverlapFringe::updateIteration is not going to be
        // called again, i.e. final overlap fringe assembly is found.

        // Iterator to iterationDataHistory_ beginning
        FIFOStack<iterationData>::iterator begin =
            iterationDataHistory_.begin();

        // Iterator to iterationDataHistory_ end
        FIFOStack<iterationData>::iterator end =
            iterationDataHistory_.end();

        // Get an iterator to iterationData object with maximum average
        // donor/acceptor suitability
        FIFOStack<iterationData>::iterator maxObject = std::max_element
        (
            begin,
            end,
            sortBySuitabilityValue
        );

        // Delete demand driven data before final holes list construction
        deleteDemandDrivenData(fringeHolesPtr_);

        // Get a reference to donor/acceptor pairs
        const donorAcceptorList& unfilteredDAPairs =
            maxObject().donorAcceptorPairs();

        // Filter donor/acceptor pairs

        // We need to clean up a bit. Namely, it is possible that a certain
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

        // Get fringe holes
        dynamicLabelList fringeHoles(maxObject().fringeHoles());

        // Loop through all fringe holes and mark them
        forAll (fringeHoles, hcI)
        {
            holeIndicatorIn[fringeHoles[hcI]] = 1.0;
        }

        // Loop through all acceptors and mark them
        forAll (unfilteredDAPairs, daPairI)
        {
            holeIndicatorIn[unfilteredDAPairs[daPairI].acceptorCell()] = 1.0;
        }

        // Get boundary field
        volScalarField::GeometricBoundaryField& holeIndicatorb =
            holeIndicator.boundaryField();

        // Perform update accross coupled boundaries, excluding overset patch
        overlapFringe::evaluateNonOversetBoundaries(holeIndicatorb);

        // Get necessary mesh data
        const cellList& meshCells = mesh.cells();
        const unallocLabelList& own = mesh.owner();
        const unallocLabelList& nei = mesh.neighbour();

        // List of acceptors to be converted to holes
        boolList accBecomingHoles(unfilteredDAPairs.size(), false);

        // Loop through all donor/acceptor pairs collected so far
        forAll (unfilteredDAPairs, daPairI)
        {
            // Get acceptor cell index
            const label& accI = unfilteredDAPairs[daPairI].acceptorCell();

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
        // converted to holes into fringeHoles and insert all other acceptors
        // into finalDAPairs temporary container
        // Create another dynamic list to collect final donor/acceptor pairs
        donorAcceptorDynamicList finalDAPairs
        (
            unfilteredDAPairs.size()
        );

        // Count number of acceptor holes that need to be converted to holes
        label nAccToHoles = 0;

        // Loop all current donor/acceptor pairs
        forAll(unfilteredDAPairs, daPairI)
        {
            if (accBecomingHoles[daPairI])
            {
                // Append the acceptor to list of holes
                fringeHoles.append(unfilteredDAPairs[daPairI].acceptorCell());
                ++nAccToHoles;
            }
            else
            {
                // Append the donor/acceptor pair to finalDAPairs list
                finalDAPairs.append(unfilteredDAPairs[daPairI]);
            }
        }

        // Print how many acceptors have been converted to holes
        reduce(nAccToHoles, sumOp<label>());
        Info<< "Converted " << nAccToHoles << " acceptors to holes." << endl;

        // Bugfix: Although we have found suitable overlap, we need to update
        // acceptors as well because eligible donors for acceptors of other
        // regions are calculated based on these acceptors (and holes)
        labelList& acceptors = *acceptorsPtr_;
        acceptors.setSize(finalDAPairs.size());
        forAll (acceptors, aI)
        {
            acceptors[aI] = finalDAPairs[aI].acceptorCell();
        }

        // Construct final donor/acceptor list
        finalDonorAcceptorsPtr_ = new donorAcceptorList
        (
            finalDAPairs.xfer()
        );

        // Construct final fringe holes list
        fringeHolesPtr_ = new labelList
        (
            fringeHoles.xfer()
        );

        // Print information
        Info<< nl
            << "Finished assembling overlap fringe for region "
            << region().name() << "." <<  nl
            << "Chosen overlap is from " << maxObject().iteration() << "."
            << " iteration. " << nl
            << "Average donor/acceptor suitability is "
            << maxObject().suitability()*100 << "%."
            << nl << endl;

        // Final overlap is found. Clear iteration history.
        iterationDataHistory_.clear();

        // Iterative process is finished. Set the flag to true.
        updateSuitableOverlapFlag(true);
    }

    return foundSuitableOverlap();
}


const Foam::labelList& Foam::adaptiveOverlapFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    // Debug, write fringe holes as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor fringe holes into a cell set for region"
            << region().name() << endl;

        cellSet holesSet
        (
            mesh(),
            "fringeHolesProc" + name(Pstream::myProcNo()) + region().name(),
            labelHashSet(*fringeHolesPtr_)
        );

        holesSet.write();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::adaptiveOverlapFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    // Debug, write candidate acceptors as a cell set
    if (oversetMesh::debug())
    {
        Pout<< "Writing processor candidate acceptors into a cell set "
            << "for region" << region().name() << endl;

        cellSet candidateAcceptorsSet
        (
            mesh(),
            "candidateAcceptorsProc" + name(Pstream::myProcNo())
          + region().name(),
            labelHashSet(*acceptorsPtr_)
        );

        candidateAcceptorsSet.write();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::adaptiveOverlapFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("adaptiveOverlapFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called adaptiveOverlapFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("adaptiveOverlapFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::adaptiveOverlapFringe::update() const
{
    Info<< "adaptiveOverlapFringe::update() const" << endl;

    // Clear everything, including acceptors and fringe holes
    clearAddressing();

    // Reset iteration counter
    fringeIter_ = 0;

    // Set flag to false
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
