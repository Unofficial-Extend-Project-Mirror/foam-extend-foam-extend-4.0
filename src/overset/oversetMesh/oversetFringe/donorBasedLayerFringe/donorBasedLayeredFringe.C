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

#include "donorBasedLayeredFringe.H"
#include "oversetRegion.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(donorBasedLayeredFringe, 0);
    addToRunTimeSelectionTable
    (
        oversetFringe,
        donorBasedLayeredFringe,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::donorBasedLayeredFringe::calcAddressing() const
{
    if (fringeHolesPtr_ || acceptorsPtr_)
    {
        FatalErrorIn
        (
            "void Foam::donorBasedLayeredFringe::calcAddressing() const"
        )   << "Fringe addressing already calculated"
            << abort(FatalError);
    }

    // We need to wait for the fringe from donor region to finish the
    // assembly. Simply create empty acceptor and fringe hole lists
    acceptorsPtr_ = new labelList();
    fringeHolesPtr_ = new labelList();
}


void Foam::donorBasedLayeredFringe::clearAddressing() const
{
    deleteDemandDrivenData(fringeHolesPtr_);
    deleteDemandDrivenData(acceptorsPtr_);
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from dictionary
Foam::donorBasedLayeredFringe::donorBasedLayeredFringe
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
    nLayers_(readLabel(dict.lookup("nLayers")))
{
    // Sanity check for number of layers
    if (nLayers_ < 1)
    {
        FatalErrorIn
        (
            "donorBasedLayeredFringe::donorBasedLayeredFringe\n"
            "(\n"
            "    const fvMesh& mesh,\n"
            "    const oversetRegion& region,\n"
            "    const dictionary& dict,\n"
            ")\n"
        )   << "Invalid number of layers specified in donor based layered"
            << "fringe dictionary."
            << nl
            << "Please specify value between greater than 0."
            << endl;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::donorBasedLayeredFringe::~donorBasedLayeredFringe()
{
    clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::donorBasedLayeredFringe::updateIteration
(
    donorAcceptorList& donorAcceptorRegionData
) const
{
    if (finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("donorBasedLayeredFringe::updateIteration(donorAcceptorList&")
            << "finalDonorAcceptorPtr_ already allocated. Something went "
            << "wrong with the iteration procedure (flag was not updated)."
            << abort(FatalError);
    }

    // Get overset mesh
    const oversetMesh& om = region().overset();

    // Get list of donor regions of this overset region
    const labelList& donorRegions = region().donorRegions();

    // Check whether all donor regions have finished with the donor search



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


const Foam::labelList& Foam::donorBasedLayeredFringe::fringeHoles() const
{
    if (!fringeHolesPtr_)
    {
        calcAddressing();
    }

    return *fringeHolesPtr_;
}


const Foam::labelList& Foam::donorBasedLayeredFringe::candidateAcceptors() const
{
    if (!acceptorsPtr_)
    {
        calcAddressing();
    }

    return *acceptorsPtr_;
}


Foam::donorAcceptorList& Foam::donorBasedLayeredFringe::finalDonorAcceptors() const
{
    if (!finalDonorAcceptorsPtr_)
    {
        FatalErrorIn("donorBasedLayeredFringe::finalDonorAcceptors()")
            << "finalDonorAcceptorPtr_ not allocated. Make sure you have "
            << "called donorBasedLayeredFringe::updateIteration() before asking for "
            << "final set of donor/acceptor pairs."
            << abort(FatalError);
    }

    if (!foundSuitableOverlap())
    {
        FatalErrorIn("donorBasedLayeredFringe::finalDonorAcceptors()")
            << "Attemted to access finalDonorAcceptors but suitable overlap "
            << "has not been found. This is not allowed. "
            << abort(FatalError);
    }

    return *finalDonorAcceptorsPtr_;
}


void Foam::donorBasedLayeredFringe::update() const
{
    if (updateFringe_)
    {
        Info<< "donorBasedLayeredFringe::update() const" << endl;

        // Clear out
        clearAddressing();
    }

    // Set flag to false and clear final donor/acceptors only
    deleteDemandDrivenData(finalDonorAcceptorsPtr_);
    updateSuitableOverlapFlag(false);
}


// ************************************************************************* //
