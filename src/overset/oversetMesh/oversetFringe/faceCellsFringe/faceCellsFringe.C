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

#include "faceCellsFringe.H"
#include "oversetRegion.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(faceCellsFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, faceCellsFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::faceCellsFringe::calcAddressing() const
{
    if (acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::faceCellsFringe::calcAddressing() const"
        )   << "Addressing already calculated"
            << abort(FatalError);
    }

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Make a hash set to collect acceptor points
    // Note: 2 patches connecting at the corner may create a duplicate,
    // which is filtered on insertion
    labelHashSet acceptorSet;

    // Find patches and mark cells
    forAll (patchNames_, nameI)
    {
        const polyPatchID curFringePatch
        (
            patchNames_[nameI],
            mesh().boundaryMesh()
        );

        if (!curFringePatch.active())
        {
            FatalErrorIn
            (
                "void faceCellsFringe::calcAddressing() const"
            )   << "Fringe patch " << patchNames_[nameI]
                << " cannot be found"
                << abort(FatalError);
        }

        const unallocLabelList& curFaceCells =
            mesh().boundaryMesh()[curFringePatch.index()].faceCells();

        forAll (curFaceCells, fcI)
        {
            // Check if cell is in region zone
            if (rcz.whichCell(curFaceCells[fcI]) > -1)
            {
                // Found acceptor
                acceptorSet.insert(curFaceCells[fcI]);
            }
        }
    }

    // Collect acceptors
    acceptorsPtr_ = new labelList(acceptorSet.sortedToc());

    // Holes are empty for this fringe
    fringeHolesPtr_ = new labelList();
}


void Foam::faceCellsFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::faceCellsFringe::faceCellsFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    patchNames_(dict.lookup("patches")),
    fringeHolesPtr_(NULL),
    acceptorsPtr_(NULL),
    finalDonorAcceptorsPtr_(NULL),
    updateFringe_
    (
        dict.lookupOrDefault<Switch>("updateAcceptors", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::faceCellsFringe::~faceCellsFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::faceCellsFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // If the donorAcceptor list has been allocated, something went wrong with
    // the iteration procedure (not-updated flag): this function has been called
    // more than once, which should not happen for faceCellsFringe
    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("faceCellsFringe::updateIteration(donorAcceptorList&")
            << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << nl << "This should not happen for faceCellsFringe."
            << abort(FatalError);
    }

    // Allocate the list by reusing the argument list
    finalDonorAcceptorsPtr_ = new donorAcceptorList
    (
        donorAcceptorRegionData,
        true
    );

    // Set the flag to true and return
    updateSuitableOverlapFlag(true);

    return foundSuitableOverlap();
}


const Foam::labelList& Foam::faceCellsFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::faceCellsFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::faceCellsFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("faceCellsFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called faceCellsFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("faceCellsFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::faceCellsFringe::update() const
{
    if (updateFringe_)
    {
        Info<< "faceCellsFringe::update() const" << endl;

        // Clear out
        clearAddressing();
    }

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
