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
    weightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::averageValueInterpolation::~averageValueInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ListScalarFieldField&
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
