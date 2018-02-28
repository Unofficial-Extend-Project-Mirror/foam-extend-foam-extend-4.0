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

#include "inverseDistanceInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(inverseDistanceInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        inverseDistanceInterpolation,
        word
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::inverseDistanceInterpolation::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "void inverseDistanceInterpolation::calcWeights() const"
        )   << "Weights already calculated."
            << abort(FatalError);
    }

    // Allocate necessary storage
    weightsPtr_ = new ListScalarFieldField(overset().regions().size());
    ListScalarFieldField& weights = *weightsPtr_;

    // Loop through all overset regions
    forAll (overset().regions(), regionI)
    {
        // Get acceptors for this region
        const donorAcceptorList& curAcceptors =
            overset().regions()[regionI].acceptors();

        // Get weights for this region
        ScalarFieldField& regionWeights = weights[regionI];
        regionWeights.setSize(curAcceptors.size());

        // Loop through acceptors of this region
        forAll (curAcceptors, aI)
        {
            // Get current donor/acceptor pair
            const donorAcceptor& da = curAcceptors[aI];

            // Get total number of donors for this acceptor
            const label nDonors = 1 + da.extendedDonorCells().size();

            // Initialise weight field to one for donors of this acceptor
            regionWeights.set
            (
                aI,
                new scalarField
                (
                    nDonors,
                    1.0
                )
            );

            // Get weights for this acceptor
            scalarField& w = regionWeights[aI];

            // Get acceptor cell centre
            const point& accCC = da.acceptorPoint();

            // Calculate master donor weight first
            w[0] /= mag(accCC - da.donorPoint()) + SMALL;

            // Calculate neighbouring donors next
            const donorAcceptor::DynamicPointList& nbrDonorPoints =
                da.extendedDonorPoints();

            forAll (nbrDonorPoints, nbrI)
            {
                // Note nbrI + 1 for subscripting because the first entry in
                // weights corresponds to the master donor
                w[nbrI + 1] /= mag(accCC - nbrDonorPoints[nbrI]) + SMALL;
            }

            // Normalise the weight field
            w /= sum(w);

        } // End for all acceptors in this region
    } // End for all regions
}


void Foam::inverseDistanceInterpolation::clearWeights() const
{
    deleteDemandDrivenData(weightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inverseDistanceInterpolation::inverseDistanceInterpolation
(
    const oversetMesh& overset,
    const word& name
)
:
    oversetInterpolation(overset, name),
    weightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceInterpolation::~inverseDistanceInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ListScalarFieldField&
Foam::inverseDistanceInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        calcWeights();
    }

    return *weightsPtr_;
}


void Foam::inverseDistanceInterpolation::update() const
{
    Info<< "inverseDistanceInterpolation::update()" << endl;

    clearWeights();
}


// ************************************************************************* //
