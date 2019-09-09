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

#include "averageValueInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(averageValueInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        averageValueInterpolation,
        word
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::averageValueInterpolation::calcWeights() const
{
    if (weightsPtr_)
    {
        FatalErrorIn
        (
            "void averageValueInterpolation::calcWeights() const"
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
            // Get total number of donors for this acceptor
            const label nDonors =
                1 + curAcceptors[aI].extendedDonorCells().size();

            // Calculate weight field as 1 divided by number of donors
            regionWeights.set
            (
                aI,
                new scalarField
                (
                    nDonors,
                    1/scalar(nDonors)
                )
            );

        } // End for all acceptors in this region
    } // End for all regions
}


void Foam::averageValueInterpolation::clearWeights() const
{
    deleteDemandDrivenData(weightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::averageValueInterpolation::averageValueInterpolation
(
    const oversetMesh& overset,
    const word& name
)
:
    oversetInterpolation(overset, name),
    weightsPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::averageValueInterpolation::~averageValueInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::scalarFieldFieldList&
Foam::averageValueInterpolation::weights() const
{
    if (!weightsPtr_)
    {
        calcWeights();
    }

    return *weightsPtr_;
}


void Foam::averageValueInterpolation::update() const
{
    Info<< "averageValueInterpolation::update()" << endl;

    clearWeights();
}


// ************************************************************************* //
