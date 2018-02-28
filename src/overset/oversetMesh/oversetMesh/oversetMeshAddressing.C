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

#include "oversetMesh.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "polyPatchID.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::oversetMesh::calcCellClassification() const
{
    if (acceptorCellsPtr_ || donorCellsPtr_ || holeCellsPtr_)
    {
        FatalErrorIn("void oversetMesh::calcCellClassification() const")
            << "Cell classification already calculated"
            << abort(FatalError);
    }

    // Count acceptor and donor and hole cells
    label nAcceptorCells = 0;
    label nDonorCells = 0;
    label nHoleCells = 0;

    forAll (regions_, regionI)
    {
        nAcceptorCells += regions_[regionI].acceptors().size();
        nDonorCells += regions_[regionI].donors().size();
        nHoleCells += regions_[regionI].holes().size();
    }

    Pout<< "Number of acceptor cells: " << nAcceptorCells << endl;
    acceptorCellsPtr_ = new labelList(nAcceptorCells);
    labelList& acceptor = *acceptorCellsPtr_;

    Pout<< "Number of donor cells: " << nDonorCells << endl;
    donorCellsPtr_ = new labelList(nDonorCells);
    labelList& donor = *donorCellsPtr_;

    Pout<< "Number of hole cells: " << nHoleCells << endl;
    holeCellsPtr_ = new labelList(nHoleCells);
    labelList& hole = *holeCellsPtr_;

    // Reset counters
    nAcceptorCells = 0;
    nDonorCells = 0;
    nHoleCells = 0;

    forAll (regions_, regionI)
    {
        // Acceptors
        const donorAcceptorList& curAcceptors = regions_[regionI].acceptors();

        forAll (curAcceptors, aI)
        {
            acceptor[nAcceptorCells] = curAcceptors[aI].acceptorCell();
            nAcceptorCells++;
        }

        // Donors
        const donorAcceptorList& curDonors = regions_[regionI].donors();

        forAll (curDonors, dI)
        {
            donor[nDonorCells] = curDonors[dI].donorCell();
            nDonorCells++;
        }

        // Holes
        const labelList& curHoles = regions_[regionI].holes();

        forAll (curHoles, hI)
        {
            hole[nHoleCells] = curHoles[hI];
            nHoleCells++;
        }
    }

    // Check donor and acceptor assembly
    Info<< "Checking donor-acceptor assembly" << endl;

    // Check for multiple acceptors

    // Check for donors that are holes

    // Check for donors that are acceptors
}


void Foam::oversetMesh::calcDomainMarkup() const
{
    if (oversetTypesPtr_ || regionIDPtr_)
    {
        FatalErrorIn("void oversetMesh::calcDomainMarkup() const")
            << "Domain markup already calculated"
            << abort(FatalError);
    }

    // Overset types
    oversetTypesPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetTypes",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimless, ACTIVE)
    );
    volScalarField& overTypes = *oversetTypesPtr_;

    scalarField& overTypesIn = overTypes.internalField();

    // Mark donor cells
    const labelList& d = donorCells();

    forAll (d, cellI)
    {
        overTypesIn[d[cellI]] = DONOR;
    }

    // Mark acceptor cells
    const labelList& a = acceptorCells();

    forAll (a, cellI)
    {
        overTypesIn[a[cellI]] = ACCEPTOR;
    }

    // Mark all hole cells
    const labelList& h = holeCells();

    forAll (h, cellI)
    {
        overTypesIn[h[cellI]] = HOLE;
    }


    // Update boundary values
    forAll (overTypes.boundaryField(), patchI)
    {
        if (!overTypes.boundaryField()[patchI].coupled())
        {
            overTypes.boundaryField()[patchI] =
                overTypes.boundaryField()[patchI].patchInternalField();
        }
    }

    // Region ID
    regionIDPtr_ = new labelList(mesh().nCells(), -1);
    labelList& rID = *regionIDPtr_;

    // Mark regions

    forAll (regions_, regionI)
    {
        const labelList& curCells = regions_[regionI].zone();

        forAll (curCells, curCellI)
        {
            rID[curCells[curCellI]] = regionI;
        }
    }

    // Check regions
    if (min(rID) < 0)
    {
        FatalErrorIn("void oversetMesh::calcDomainMarkup() const")
            << "Found cells without region ID.  Please check overset setup"
            << abort(FatalError);
    }
}


void Foam::oversetMesh::calcGamma() const
{
    if (gammaPtr_ || gammaExtPtr_ || sGammaPtr_)
    {
        FatalErrorIn("void oversetMesh::calcGamma() const")
            << "Markup fields already calculated"
            << abort(FatalError);
    }

    // Trigger the calculation of cell centres to avoid tangled parallel
    // communications in case the lazy evaluation mechanism is invoked in
    // GeometricField::evaluate() member function. Temporary solution.
    // To-do: examine stack trace for fvMesh::makeC(), VV, 8/Feb/2016.
    mesh().C();

    // Fluid cells indicator, marking only live cells
    gammaPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetGamma",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1)
    );
    volScalarField& g = *gammaPtr_;
    scalarField& gIn = g.internalField();

    // Non-hole cells indicator, marking live and acceptor cells
    gammaExtPtr_ = new volScalarField
    (
        IOobject
        (
            "oversetGammaExt",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 1)
    );
    volScalarField& gExt = *gammaExtPtr_;
    scalarField& gExtIn = gExt.internalField();

    // Mark hole cells
    const labelList& holes = holeCells();

    forAll (holes, holeI)
    {
        gIn[holes[holeI]] = 0;
        gExtIn[holes[holeI]] = 0;
    }

    // Mark acceptor cells: they are marked with 1 in gammaExt, but with
    // 0 in gamma
    const labelList& acceptors = acceptorCells();

    forAll (acceptors, accI)
    {
        gIn[acceptors[accI]] = 0;
    }

    // Before collecting active faces, update coupled boundaries
    // to make sure patchNeighbourField is correct

    // Not allowed to call correctBoundaryConditions.  HJ, 16/Apr/2012
    // Evaluate coupled boundaries and copy out the uncoupled ones

    g.boundaryField().evaluateCoupled();

    forAll (g.boundaryField(), patchI)
    {
        if (!g.boundaryField()[patchI].coupled())
        {
            g.boundaryField()[patchI] =
                g.boundaryField()[patchI].patchInternalField();
        }
    }

    gExt.boundaryField().evaluateCoupled();

    forAll (gExt.boundaryField(), patchI)
    {
        if (!gExt.boundaryField()[patchI].coupled())
        {
            gExt.boundaryField()[patchI] =
                gExt.boundaryField()[patchI].patchInternalField();
        }
    }

    // Calculate face mask
    sGammaPtr_ = new surfaceScalarField
    (
        IOobject
        (
            "oversetSGamma",
            mesh().time().timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("one", dimless, 0)
    );
    surfaceScalarField& sg = *sGammaPtr_;
    scalarField& sgIn = sg.internalField();

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    // Internal faces: flux is live between all active and acceptor cells
    forAll (sgIn, faceI)
    {
        // If both cells are live, flux is live
        if
        (
            gExtIn[owner[faceI]] > SMALL
         && gExtIn[neighbour[faceI]] > SMALL
        )
        {
            sgIn[faceI] = 1;
        }
    }

    // Remove faces where both owner and neighbour are acceptors
    volScalarField gAcc = gExt - g;
    gAcc.boundaryField().evaluateCoupled();

    const scalarField& gAccIn = gAcc.internalField();

    // Kill fluxes between two acceptor cells
    forAll (sgIn, faceI)
    {
        // If both cells are live, flux is live
        if
        (
            gAccIn[owner[faceI]] > SMALL
         && gAccIn[neighbour[faceI]] > SMALL
        )
        {
            sgIn[faceI] = 0;
        }
    }

    surfaceScalarField::GeometricBoundaryField& sgPatches =
        sg.boundaryField();

    const volScalarField::GeometricBoundaryField& gExtPatches =
        gExt.boundaryField();

    const volScalarField::GeometricBoundaryField& gAccPatches =
        gAcc.boundaryField();


    forAll (gExtPatches, patchI)
    {
        if (gExtPatches[patchI].coupled())
        {
            scalarField& gP = sgPatches[patchI];

            // For coupled patches, check gammaExt
            const scalarField gammaOwn =
                gExtPatches[patchI].patchInternalField();

            const scalarField gammaNei =
                gExtPatches[patchI].patchNeighbourField();

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

            const scalarField gammaAccOwn =
                gAccPatches[patchI].patchInternalField();

            const scalarField gammaAccNei =
                gAccPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, faceI)
            {
                if
                (
                    gammaAccOwn[faceI] > SMALL
                 && gammaAccNei[faceI] > SMALL
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
            const scalarField gammaFc =
                g.boundaryField()[patchI].patchInternalField();

            scalarField& gP = sgPatches[patchI];

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


void Foam::oversetMesh::calcFringeFaces() const
{
    if
    (
        fringeFacesPtr_
     || fringeFaceCellsPtr_
     || fringeFaceFlipsPtr_
     || acceptorInternalFacesPtr_
    )
    {
        FatalErrorIn("void oversetMesh::calcFringeFaces() const")
            << "Face masks already calculated"
            << abort(FatalError);
    }

    labelList acCellIndicator(mesh().nCells(), -1);

    // Mark hole cells with -2
    const labelList& hc = holeCells();

    forAll (hc, hcI)
    {
        acCellIndicator[hc[hcI]] = -2;
    }

    // Mark acceptor cells with their index
    const labelList& acc = acceptorCells();

    forAll (acc, accI)
    {
        acCellIndicator[acc[accI]] = accI;
    }

    DynamicList<label> acF(2*acc.size());
    DynamicList<label> acFC(2*acc.size());
    DynamicList<bool> acFF(2*acc.size());

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    const volScalarField::GeometricBoundaryField& gammaPatches =
        gamma().boundaryField();

    // Fringe face collection

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            acCellIndicator[owner[faceI]] == -1
         && acCellIndicator[neighbour[faceI]] > -1
        )
        {
            // Owner is live, neighbour AC.  Its AC index is in
            // acCellIndicator
            acF.append(faceI);
            acFC.append(acCellIndicator[neighbour[faceI]]);
            acFF.append(false);
        }
        else if
        (
            acCellIndicator[owner[faceI]] > -1
         && acCellIndicator[neighbour[faceI]] == -1
        )
        {
            // Neighbour is live, owner AC.  Its AC index is in
            // acCellIndicator
            acF.append(faceI);
            acFC.append(acCellIndicator[owner[faceI]]);
            acFF.append(true);
        }
    }

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh().boundary()[patchI].faceCells();
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {
            const scalarField gammaOwn =
                gammaPatches[patchI].patchInternalField();

            const scalarField gammaNei =
                gammaPatches[patchI].patchNeighbourField();

            forAll (gammaOwn, patchFaceI)
            {
                // Note: we assume that there are no live cells neighbouring
                // hole cells, only acceptors. This is forced in overset fringe
                // algorithms (overlapFringe specifically).
                if
                (
                    mag(gammaNei[patchFaceI] - gammaOwn[patchFaceI]) > SMALL
                )
                {
                    if (acCellIndicator[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is AC
                        acF.append(start + patchFaceI);
                        acFC.append(acCellIndicator[fc[patchFaceI]]);
                        acFF.append(true);
                    }
                    else
                    {
                        // Neighbour cell is AC
                        acF.append(start + patchFaceI);
                        acFC.append(-1);
                        acFF.append(false);
                    }
                }
            }
        }
    }

    // Pack the data
    acF.shrink();
    acFC.shrink();
    acFF.shrink();

    fringeFacesPtr_ = new labelList(acF.xfer());

    fringeFaceCellsPtr_ = new labelList(acFC.xfer());
    fringeFaceFlipsPtr_ = new boolList(acFF.xfer());

    // Acceptor internal face collection

    DynamicList<label> acInternalF(2*acc.size());

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            acCellIndicator[owner[faceI]] > -1
         && acCellIndicator[neighbour[faceI]] > -1
        )
        {
            // Found face between two acceptor cells
            acInternalF.append(faceI);
        }
    }

    const volScalarField::GeometricBoundaryField& gammaExtPatches =
        gammaExt().boundaryField();

    forAll (gammaPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaPatches[patchI].coupled())
        {
            // Create acceptor mask
            const scalarField gammaAccOwn =
                gammaExtPatches[patchI].patchInternalField()
              - gammaPatches[patchI].patchInternalField();

            const scalarField gammaAccNei =
                gammaExtPatches[patchI].patchNeighbourField()
              - gammaPatches[patchI].patchNeighbourField();

            forAll (gammaAccOwn, patchFaceI)
            {
                if
                (
                    mag(gammaAccNei[patchFaceI]) > SMALL
                 && mag(gammaAccOwn[patchFaceI]) > SMALL
                )
                {
                    acInternalF.append(start + patchFaceI);
                }
            }
        }
    }

    // Pack the data
    acInternalF.shrink();

    acceptorInternalFacesPtr_ = new labelList(acInternalF.xfer());
}


void Foam::oversetMesh::calcHoleFaces() const
{
    if
    (
        holeFacesPtr_
     || holeFaceCellsPtr_
     || holeFaceFlipsPtr_
     || holeInternalFacesPtr_
    )
    {
        FatalErrorIn("void oversetMesh::calcHoleFaces() const")
            << "Face masks already calculated"
            << abort(FatalError);
    }

    labelList hole(mesh().nCells(), -1);

    // Mark hole cells with their index
    const labelList& hc = holeCells();

    forAll (hc, hcI)
    {
        hole[hc[hcI]] = hcI;
    }

    DynamicList<label> hF(2*hc.size());
    DynamicList<label> hFC(2*hc.size());
    DynamicList<bool> hFF(2*hc.size());

    const unallocLabelList& owner = mesh().owner();
    const unallocLabelList& neighbour = mesh().neighbour();

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            hole[owner[faceI]] == -1
         && hole[neighbour[faceI]] > -1
        )
        {
            // Owner is live (acceptor), neighbour H.  Its H index is in
            // hole
            hF.append(faceI);
            hFC.append(hole[neighbour[faceI]]);
            hFF.append(false);
        }
        else if
        (
            hole[owner[faceI]] > -1
         && hole[neighbour[faceI]] == -1
        )
        {
            // Neighbour is live (acceptor), owner H.  Its H index is in
            // hole
            hF.append(faceI);
            hFC.append(hole[owner[faceI]]);
            hFF.append(true);
        }
    }

    // Bug fix: need to use gammaExt instead of gamma to mark hole faces
    const volScalarField::GeometricBoundaryField& gammaExtPatches =
        gammaExt().boundaryField();

    forAll (gammaExtPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const labelList& fc = mesh().boundary()[patchI].faceCells();
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaExtPatches[patchI].coupled())
        {
            const scalarField gammaExtOwn =
                gammaExtPatches[patchI].patchInternalField();

            const scalarField gammaExtNei =
                gammaExtPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                // Note: we assume that there are no live cells neighbouring
                // hole cells, only acceptors. This is forced in overset fringe
                // algorithms (overlapFringe specifically).
                if
                (
                    mag(gammaExtNei[patchFaceI] - gammaExtOwn[patchFaceI])
                  > SMALL
                )
                {
                    if (hole[fc[patchFaceI]] > -1)
                    {
                        // Owner cell is H
                        hF.append(start + patchFaceI);
                        hFC.append(hole[fc[patchFaceI]]);
                        hFF.append(true);
                    }
                    else
                    {
                        // Neighbour cell is H
                        hF.append(start + patchFaceI);
                        hFC.append(-1);
                        hFF.append(false);
                    }
                }
            }
        }
    }

    // Pack the data
    hF.shrink();
    hFC.shrink();
    hFF.shrink();

    holeFacesPtr_ = new labelList(hF.xfer());

    holeFaceCellsPtr_ = new labelList(hFC.xfer());
    holeFaceFlipsPtr_ = new boolList(hFF.xfer());


    // Hole internal faces

    // Hole faces are the ones where both owner and neighbour is a hole

    DynamicList<label> hIF(2*hc.size());

    forAll (neighbour, faceI)
    {
        // Old check done using mag gamma.  Changed for internal cells
        if
        (
            hole[owner[faceI]] > -1
         && hole[neighbour[faceI]] > -1
        )
        {
            hIF.append(faceI);
        }
    }

    forAll (gammaExtPatches, patchI)
    {
        // Note: take faceCells from fvPatch (because of empty)
        const label start = mesh().boundaryMesh()[patchI].start();

        if (gammaExtPatches[patchI].coupled())
        {
            const scalarField gammaExtOwn =
                gammaExtPatches[patchI].patchInternalField();

            const scalarField gammaExtNei =
                gammaExtPatches[patchI].patchNeighbourField();

            forAll (gammaExtOwn, patchFaceI)
            {
                if
                (
                    gammaExtNei[patchFaceI] < SMALL
                 && gammaExtOwn[patchFaceI] < SMALL
                )
                {
                    hIF.append(start + patchFaceI);
                }
            }
        }
    }

    // Pack the data
    hIF.shrink();
    holeInternalFacesPtr_ = new labelList(hIF.xfer());
}


void Foam::oversetMesh::calcInterpolationMap() const
{
    if (mapPtr_ || remoteDonorToLocalAcceptorAddrPtr_)
    {
        FatalErrorIn("void oversetMesh::calcInterpolationMap() const")
            << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // Create list containing number of donors my processor is sending to other
    // processors.
    // Example: nDonorsToProcessorMap[2][5] = 123 will (after collection of
    // data) mean that processor 2 sends 123 donor values to processor 5.
    // Note that for a serial run, this list is completely
    // unnecessary, but I prefer writing this in a general way, where I don't
    // care about minor loss of efficiency for serial runs. VV, 21/Oct/2016.
    labelListList nDonorsToProcessorMap(Pstream::nProcs());

    forAll (nDonorsToProcessorMap, procI)
    {
        nDonorsToProcessorMap[procI].setSize(Pstream::nProcs(), 0);
    }

    // Algorithm:
    // 1) First, we need to calculate the sending map, which tells me which
    //    donor values I need to send to which processor:
    //    - Loop through all regions and then through all donors for that region
    //    - Mark all donors associated with donor/acceptor pair and see to which
    //      processor I need to send its data. Handle the possiblity that a
    //      single donor may be used by multiple acceptors on the same processor
    //    - While looping through donors, count how many donors I'm sending to
    //      each processor
    //    - While I'm calculating the sending map, I will keep track of mapping
    //      from local donor cell indices to actually sent indices. This will
    //      allow me to easily create remoteDonorToLocalAcceptorAddressing
    // 2) Second, we need to count how many items my (local) processor is going
    //    to be receiving from all other processors
    // 3) Third, I need to create constructing map, which will be organized as
    //    follows:
    //    - If processor N sends me M donor values, these values will be stored
    //      in the receiving list from indices K to K + M - 1, where K is the
    //      total number of donors received from previous processors (processors
    //      I, where I < N).

    // STAGE 1: calculation of sending map and receiving sizes

    // Get the nDonorsToProcessorMap for my processor
    labelList& numberOfLocalDonorsToProcs =
        nDonorsToProcessorMap[Pstream::myProcNo()];

    // Create sending map and roughly initialize the capacity as 6 times the
    // number of donors for this processor (based on the assumption that all
    // data goes to a single processor and we have hexahedral cells). This
    // assumption should be fine considered the memory penalty of the CFD run
    // compared to overset assembly.
    // Note: hash set is used to send only unique donors as single donor can be
    // used for multiple acceptors.
    List<labelHashSet> sendMap(Pstream::nProcs());

    // Create labelList containg donor indices which will be used to create
    // remote donor to local acceptor addressing
    labelList donorIDs(mesh().nCells(), -1);

    // Enclosed in scope for readability
    {
        // Count number of local donors
        label nDonorsToSend = 0;

        forAll (regions_, regionI)
        {
            nDonorsToSend += regions_[regionI].donors().size();
        }
        nDonorsToSend *= 6;

        // Set the capacity of each sendMap (for each processor)
        forAll (sendMap, receivingProcI)
        {
            sendMap[receivingProcI].resize(nDonorsToSend);
        }
    }

    // Loop through all regions
    forAll (regions_, regionI)
    {
        // Get the list of donor/acceptor data
        const donorAcceptorList& regionDonors = regions_[regionI].donors();

        // Loop through all donors
        forAll (regionDonors, dI)
        {
            // Get necessary master donor data
            const donorAcceptor& curDonor = regionDonors[dI];
            const label donorCellI = curDonor.donorCell();
            const label acceptorProcIndex = curDonor.acceptorProcNo();

            // Set donor index
            donorIDs[donorCellI] = donorCellI;

            // Insert this master donor into the hash table for the processor
            // I'm sending data to
            if (sendMap[acceptorProcIndex].insert(donorCellI))
            {
                // This pair has not been registered yet, increment the counter
                ++numberOfLocalDonorsToProcs[acceptorProcIndex];
            }

            // Insert extended donors for this donor/acceptor pair
            const donorAcceptor::DynamicLabelList& extDonors =
                curDonor.extendedDonorCells();

            forAll (extDonors, eDonorCellI)
            {
                // Get index
                const label extDonorI = extDonors[eDonorCellI];

                // Set donor index for this extended donor
                donorIDs[extDonorI] = extDonorI;

                if (sendMap[acceptorProcIndex].insert(extDonorI))
                {
                    ++numberOfLocalDonorsToProcs[acceptorProcIndex];
                }
            }
        } // End for all donors in this region
    } // End for all regions

    // Create sending map from the List of hash sets
    labelListList sendDataMap(Pstream::nProcs());
    forAll (sendDataMap, procI)
    {
        sendDataMap[procI] = sendMap[procI].toc();
    }

    // Gather/scatter number of donors going to each processor from each
    // processor so that all processors have all necessary information when
    // creating map distribute tool
    Pstream::gatherList(nDonorsToProcessorMap);
    Pstream::scatterList(nDonorsToProcessorMap);

    // For my processor, I need to count how many items I'm going to be
    // receiving from others. Note that the receive size is different for each
    // processor (we are not collecting data in a global list as in e.g. GGI)
    label nReceives = 0;
    forAll (nDonorsToProcessorMap, procI)
    {
        nReceives += nDonorsToProcessorMap[procI][Pstream::myProcNo()];
    }

    // STAGE 2: calculation of construct map

    // The construct map is simply an offset of the index by a number of values
    // received by previous processors
    // Example:
    /*
       Procs sending to me | Number of items being sent
      --------------------------------------------------
              P0           |             1              
              P1           |             7              
              P5           |             2              
              .            |             .
              .            |             .
              .            |             .

        Received data has the following form:
        (
            a_0, (one value from proc 0)
            a_1, a_2, a_3, a_4, a_5, a_6, a_7 (seven values from proc 1)
            a_8, a_9 (two values from proc 5)
            ...
            ...
            ...
    */

    // Create construct map
    labelListList constructDataMap(Pstream::nProcs());

    // Counter for offset
    label procOffset = 0;

    forAll (constructDataMap, procI)
    {
        // Get receiving size from this processor
        const label nReceivesFromCurProc =
            nDonorsToProcessorMap[procI][Pstream::myProcNo()];

        // Get current construct map
        labelList& curConstructMap = constructDataMap[procI];

        // Set the size of each list corresponding to number of received values
        // from this processor
        curConstructMap.setSize(nReceivesFromCurProc);

        // Set mapping via simple offset
        forAll (curConstructMap, receivedItemI)
        {
            curConstructMap[receivedItemI] = receivedItemI + procOffset;
        }

        // Increment the processor offset by the size received from this
        // processor
        procOffset += nReceivesFromCurProc;
    }

    // STAGE 3: create the map distribute object

    // Create the map distribute object
    mapPtr_ = new mapDistribute
    (
        nReceives,
        sendDataMap,
        constructDataMap
    );

    // STAGE 4: create remoteDonorToLocalAcceptor addressing object

    // Distribute local donor cell ID's using mapDistribute
    mapPtr_->distribute(donorIDs);

    // Sanity check whether all donorIDs that have been sent to me are valid
    if (min(donorIDs) < 0)
    {
        FatalErrorIn("void oversetMesh::calcInterpolationMap() const")
            << "Found invalid donor index after distribution." << nl
            << "Something has seriously gone wrong..."
            << abort(FatalError);
    }

    // Now we have all received donorIDs from all processors, combined in a
    // single list. We know which part has come from which processor, and we
    // need to disect this data in order to create a list (indexed by processor
    // ID) of labelFields (indexed by cell index), which will eventually give us
    // the index of donor in the received list (for a given donor processor ID
    // and local donor index)
    remoteDonorToLocalAcceptorAddrPtr_ =
        new List<labelField>(Pstream::nProcs());
    List<labelField>& rdlaAddr = *remoteDonorToLocalAcceptorAddrPtr_;

    // Each processor needs to know how many cells we have in all other
    // processors
    labelList procNCells(Pstream::nProcs());
    procNCells[Pstream::myProcNo()] = mesh().nCells();

    Pstream::gatherList(procNCells);
    Pstream::scatterList(procNCells);

    // Allocate storage and initialise to -1. Note that each processor stores
    // nProcs fields of the size corresponding to local processor mesh size
    // (nCells). This is a bit wasteful, but allows us efficient interpolation
    forAll (rdlaAddr, procI)
    {
        rdlaAddr[procI].setSize(procNCells[procI], -1);
    }

    // Helper: starting processor index for list slicing
    label startProcIndex = 0;

    // Loop through processors
    for (label procI = 0; procI < Pstream::nProcs(); ++procI)
    {
        // Get number of entries I have received from this processor
        const label nCurProcRec =
            nDonorsToProcessorMap[procI][Pstream::myProcNo()];

        // Get the corresponding subList
        const labelList::subList procDonorIDs
        (
            donorIDs, // Original list
            nCurProcRec, // Size of the data
            startProcIndex // Start index
        );

        // Fill in necessary parts of the remote donor to local acceptor
        // addressing
        labelField& curAddr = rdlaAddr[procI];
        forAll (procDonorIDs, i)
        {
            // Offset the index
            label offsetI = i + startProcIndex;

            // Store distributed donorID in the location
            curAddr[procDonorIDs[i]] = offsetI;
        }

        // Increment the startProcIndex by number of receives from this
        // processor
        startProcIndex += nCurProcRec;
    }


    // Force calculation of domain markup fields for post-processing
    // HJ, 9/Apr/2013
    oversetTypes();
}


void Foam::oversetMesh::clearOut() const
{
    deleteDemandDrivenData(acceptorCellsPtr_);
    deleteDemandDrivenData(donorCellsPtr_);
    deleteDemandDrivenData(holeCellsPtr_);

    deleteDemandDrivenData(oversetTypesPtr_);
    deleteDemandDrivenData(regionIDPtr_);

    deleteDemandDrivenData(gammaPtr_);
    deleteDemandDrivenData(gammaExtPtr_);
    deleteDemandDrivenData(sGammaPtr_);

    deleteDemandDrivenData(fringeFacesPtr_);
    deleteDemandDrivenData(fringeFaceCellsPtr_);
    deleteDemandDrivenData(fringeFaceFlipsPtr_);

    deleteDemandDrivenData(holeFacesPtr_);
    deleteDemandDrivenData(holeFaceCellsPtr_);
    deleteDemandDrivenData(holeFaceFlipsPtr_);
    deleteDemandDrivenData(holeInternalFacesPtr_);
    deleteDemandDrivenData(acceptorInternalFacesPtr_);

    deleteDemandDrivenData(mapPtr_);
    deleteDemandDrivenData(remoteDonorToLocalAcceptorAddrPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelList& Foam::oversetMesh::acceptorCells() const
{
    if (!acceptorCellsPtr_)
    {
        calcCellClassification();
    }

    return *acceptorCellsPtr_;
}


const Foam::labelList& Foam::oversetMesh::donorCells() const
{
    if (!donorCellsPtr_)
    {
        calcCellClassification();
    }

    return *donorCellsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeCells() const
{
    if (!holeCellsPtr_)
    {
        calcCellClassification();
    }

    return *holeCellsPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::oversetTypes() const
{
    if (!oversetTypesPtr_)
    {
        calcDomainMarkup();
    }

    return *oversetTypesPtr_;
}


const Foam::labelList& Foam::oversetMesh::regionID() const
{
    if (!regionIDPtr_)
    {
        calcDomainMarkup();
    }

    return *regionIDPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::gamma() const
{
    if (!gammaPtr_)
    {
        calcGamma();
    }

    return *gammaPtr_;
}


const Foam::volScalarField& Foam::oversetMesh::gammaExt() const
{
    if (!gammaExtPtr_)
    {
        calcGamma();
    }

    return *gammaExtPtr_;
}


const Foam::surfaceScalarField& Foam::oversetMesh::sGamma() const
{
    if (!sGammaPtr_)
    {
        calcGamma();
    }

    return *sGammaPtr_;
}


const Foam::labelList& Foam::oversetMesh::fringeFaces() const
{
    if (!fringeFacesPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::fringeFaceCells() const
{
    if (!fringeFaceCellsPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFaceCellsPtr_;
}


const Foam::boolList& Foam::oversetMesh::fringeFaceFlips() const
{
    if (!fringeFaceFlipsPtr_)
    {
        calcFringeFaces();
    }

    return *fringeFaceFlipsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeFaces() const
{
    if (!holeFacesPtr_)
    {
        calcHoleFaces();
    }

    return *holeFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeFaceCells() const
{
    if (!holeFaceCellsPtr_)
    {
        calcHoleFaces();
    }

    return *holeFaceCellsPtr_;
}


const Foam::boolList& Foam::oversetMesh::holeFaceFlips() const
{
    if (!holeFaceFlipsPtr_)
    {
        calcHoleFaces();
    }

    return *holeFaceFlipsPtr_;
}


const Foam::labelList& Foam::oversetMesh::holeInternalFaces() const
{
    if (!holeInternalFacesPtr_)
    {
        calcHoleFaces();
    }

    return *holeInternalFacesPtr_;
}


const Foam::labelList& Foam::oversetMesh::acceptorInternalFaces() const
{
    if (!acceptorInternalFacesPtr_)
    {
        calcFringeFaces();
    }

    return *acceptorInternalFacesPtr_;
}


const Foam::mapDistribute& Foam::oversetMesh::map() const
{
    if (!mapPtr_)
    {
        calcInterpolationMap();
    }

    return *mapPtr_;
}


const Foam::List<Foam::labelField>&
Foam::oversetMesh::remoteDonorToLocalAcceptorAddr() const
{
    if (!remoteDonorToLocalAcceptorAddrPtr_)
    {
        calcInterpolationMap();
    }

    return *remoteDonorToLocalAcceptorAddrPtr_;
}


// ************************************************************************* //
