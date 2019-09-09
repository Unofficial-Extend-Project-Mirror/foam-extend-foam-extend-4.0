/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    weightsPtr_ = new scalarFieldFieldList(overset().regions().size());
    scalarFieldFieldList& weights = *weightsPtr_;

    // Loop through all overset regions
    forAll (overset().regions(), regionI)
    {
        // Get acceptors for this region
        const donorAcceptorList& curAcceptors =
            overset().regions()[regionI].acceptors();

        // Get weights for this region
        scalarFieldField& regionWeights = weights[regionI];
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
    weightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inverseDistanceInterpolation::~inverseDistanceInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::scalarFieldFieldList&
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
