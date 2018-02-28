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

#include "manualFringe.H"
#include "oversetRegion.H"
#include "foamTime.H"
#include "cellSet.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(manualFringe, 0);
    addToRunTimeSelectionTable(oversetFringe, manualFringe, dictionary);
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::manualFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::manualFringe::calcAddressing() const"
        )   << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    fringeHolesPtr_ = new labelList
    (
        cellSet
        (
            mesh(),
            holesSetName_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ).toc()
    );

    acceptorsPtr_ = new labelList
    (
        cellSet
        (
            mesh(),
            acceptorsSetName_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ).toc()
    );

    // Debug

    // Get reference to region cell zone
    const cellZone& rcz = region().zone();

    // Check holes
    const labelList& h = *fringeHolesPtr_;

    forAll (h, holeI)
    {
        if (rcz.whichCell(h[holeI]) < 0)
        {
            FatalErrorIn
            (
                "void Foam::manualFringe::calcAddressing() const"
            )   << "Invalid hole cell for region " << region().name()
                << ": cell " << h[holeI] << " does not belong to this region"
                << abort(FatalError);
        }
    }

    // Check acceptors
    const labelList& a = *acceptorsPtr_;

    forAll (a, accI)
    {
        if (rcz.whichCell(a[accI]) < 0)
        {
            FatalErrorIn
            (
                "void Foam::manualFringe::calcAddressing() const"
            )   << "Invalid acceptor cell for region " << region().name()
                << ": cell " << a[accI] << " does not belong to this region"
                << abort(FatalError);
        }
    }
}


void Foam::manualFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::manualFringe::manualFringe
(
    const fvMesh& mesh,
    const oversetRegion& region,
    const dictionary& dict
)
:
    oversetFringe(mesh, region, dict),
    holesSetName_(dict.lookup("holes")),
    acceptorsSetName_(dict.lookup("acceptors")),
    fringeHolesPtr_(NULL),
    acceptorsPtr_(NULL),
    finalDonorAcceptorsPtr_(NULL),
    updateFringe_
    (
        dict.lookupOrDefault<Switch>("updateAcceptors", false)
    )
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::manualFringe::~manualFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::manualFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    // If the donorAcceptor list has been allocated, something went wrong with
    // the iteration procedure (not-updated flag): this function has been called
    // more than once, which should not happen for manualFringe
    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("manualFringe::updateIteration(donorAcceptorList&)")
            << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << nl << "This should not happen for manualFringe."
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


const Foam::labelList& Foam::manualFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::manualFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::manualFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("manualFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called manualFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("manualFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::manualFringe::update() const
{
    if (updateFringe_)
    {
        Info<< "manualFringe::update() const" << endl;

        // Clear out
        clearAddressing();
    }

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
