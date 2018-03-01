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

#include "leastSquaresGradientInterpolation.H"
#include "oversetInterpolation.H"
#include "oversetMesh.H"
#include "leastSquaresVectors.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(leastSquaresGradientInterpolation, 0);
    addToRunTimeSelectionTable
    (
        oversetInterpolation,
        leastSquaresGradientInterpolation,
        dictionary
    );
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::leastSquaresGradientInterpolation::calcWeights() const
{
    if (localWeightsPtr_ || remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "void leastSquaresGradientInterpolation::calcWeights() const"
        )   << "Weights already calculated."
            << abort(FatalError);
    }

    // Get list of local donors
    const labelList& ld = overset().localDonors();

    // Get list of neighbouring donors (needed to set appropriate size of each
    // bottom most list of weights - for all donors)
    const labelListList& lnd = overset().localNeighbouringDonors();

    // Get local donor addressing (local acceptor index for each local donor)
    const labelList& ldAddr = overset().localDonorAddr();

    // Get local acceptors and donors
    const labelList& acceptorCells = overset().acceptorCells();

    // Get necessary mesh data
    const vectorField& CC = overset().mesh().cellCentres();
    const lduAddressing& lduAddr = overset().mesh().lduAddr();

    // Get least squares vectors
    const leastSquaresVectors& lsv =
        leastSquaresVectors::New(overset().mesh());
    const vectorField& ownLsIn = lsv.pVectors().internalField();
    const vectorField& neiLsIn = lsv.nVectors().internalField();

    // Create a local weight field
    localWeightsPtr_ = new ScalarFieldField(ld.size());
    ScalarFieldField& localWeights = *localWeightsPtr_;

    // Loop through local donors, setting the weights
    forAll (ld, ldI)
    {
        // Set the size of this donor weight field (master donor (1) +
        // neighbouring donors). Set the size and all values to 0.
        localWeights.set
        (
            ldI,
            new scalarField(1 + lnd[ldI].size(), 0)
        );

        // Get acceptor cell centre
        const vector& accCC = CC[acceptorCells[ldAddr[ldI]]];

        // Least squares gradient interpolation: weights are defined using
        // extrapolation from the master donor cell based on least squares
        // gradient

        // Get reference to current weight field for this acceptor
        scalarField& curLocWeights = localWeights[ldI];

        // Set master donor weight to one
        scalar& masterDonorWeight = curLocWeights[0];
        masterDonorWeight = 1;

        // Get master donor index
        const label& masterDonorI = ld[ldI];

        // Calculate distance vector from master donor to acceptor
        const vector dDA = accCC - CC[masterDonorI];

        // Calculate neighbouring donors next
        const labelList& curNbrDonors = lnd[ldI];
        forAll (curNbrDonors, nbrI)
        {
            // Get neighbouring donor index
            const label nbrDonorI = curNbrDonors[nbrI];

            // Get face index from addressing
            const label faceI = lduAddr.triIndex(masterDonorI, nbrDonorI);

            // Get reference to the weight for this neighbouring donor
            scalar& curNbrWeight = curLocWeights[nbrI + 1];

            if (masterDonorI < nbrDonorI)
            {
                // Master is owner, calculate neighbouring donor weight using p
                // vector at this face
                curNbrWeight = dDA & ownLsIn[faceI];
            }
            else
            {
                // Master is neighbour, calculate neighbouring donor weight
                // using n vector at this face
                curNbrWeight = dDA & neiLsIn[faceI];
            }

            // Update master donor weight
            masterDonorWeight -= curNbrWeight;
        }

        // Note: weight field normalized by definition of least squares
        // interpolation
    }

    // Handling remote donors for a parallel run
    if (Pstream::parRun())
    {
        // Get remote donor addressing
        const labelList& rd = overset().remoteDonors();
        const labelListList& rnd = overset().remoteNeighbouringDonors();

        // Create a global weights list for remote donor weights
        remoteWeightsPtr_ = new ListScalarFieldField(Pstream::nProcs());
        ListScalarFieldField& remoteWeights = *remoteWeightsPtr_;

        // Get the corresponding acceptor cell centres for remote donors on this
        // processor
        const vectorField& procRemAccCC = remoteAccCC()[Pstream::myProcNo()];

        // Set the size of the weight field for this processor
        ScalarFieldField& myProcRemoteWeights =
            remoteWeights[Pstream::myProcNo()];
        myProcRemoteWeights.setSize(rd.size());

        // Loop through remote donors and calculate weights
        forAll (rd, rdI)
        {
            // Allocate the storage for this donor weight field (master donor
            // (1) + neighbouring donors). Set the size and all values to 0.
            myProcRemoteWeights.set
            (
                rdI,
                new scalarField(1 + rnd[rdI].size(), 1)
            );

            // Get acceptor cell centre
            const vector& accCC = procRemAccCC[rdI];

            // Get reference to current weight field for this acceptor
            scalarField& curRemWeights = myProcRemoteWeights[rdI];

            // Set master donor weight to one
            scalar& masterDonorWeight = curRemWeights[0];
            masterDonorWeight = 1;

            // Get master donor index
            const label& masterDonorI = rd[rdI];

            // Calculate distance vector from master donor to acceptor
            const vector dDA = accCC - CC[masterDonorI];

            // Calculate neighbouring donors next
            const labelList& curNbrDonors = rnd[rdI];
            forAll (curNbrDonors, nbrI)
            {
                // Get neighbouring donor index
                const label nbrDonorI = curNbrDonors[nbrI];

                // Get face index from addressing
                const label faceI = lduAddr.triIndex(masterDonorI, nbrDonorI);

                // Get reference to the weight for this neighbouring donor
                scalar& curNbrWeight = curRemWeights[nbrI + 1];

                if (masterDonorI < nbrDonorI)
                {
                    // Master is owner, calculate neighbouring donor weight
                    // using p vector at this face
                    curNbrWeight = dDA & ownLsIn[faceI];
                }
                else
                {
                    // Master is neighbour, calculate neighbouring donor weight
                    // using n vector at this face
                    curNbrWeight = dDA & neiLsIn[faceI];
                }

                // Update master donor weight
                masterDonorWeight -= curNbrWeight;
            }

            // Note: weight field normalized by definition of least squares
            // interpolation
        }

        // Gather remote weights (no need to scatter since the data is needed
        // only for the master processor).
        Pstream::gatherList(remoteWeights);
    }
}


void Foam::leastSquaresGradientInterpolation::clearWeights() const
{
    deleteDemandDrivenData(localWeightsPtr_);
    deleteDemandDrivenData(remoteWeightsPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::leastSquaresGradientInterpolation::leastSquaresGradientInterpolation
(
    const oversetMesh& overset,
    const dictionary& dict
)
:
    oversetInterpolation(overset, dict),
    localWeightsPtr_(NULL),
    remoteWeightsPtr_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::leastSquaresGradientInterpolation::~leastSquaresGradientInterpolation()
{
    clearWeights();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::oversetInterpolation::ScalarFieldField&
Foam::leastSquaresGradientInterpolation::localWeights() const
{
    if (!localWeightsPtr_)
    {
        calcWeights();
    }

    return *localWeightsPtr_;
}


const Foam::oversetInterpolation::ListScalarFieldField&
Foam::leastSquaresGradientInterpolation::remoteWeights() const
{
    // We cannot calculate the remoteWeights using usual lazy evaluation
    // mechanism since the data only exists on the master processor. Add
    // additional guards, VV, 8/Feb/2016.
    if (!Pstream::parRun())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a serial run."
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!Pstream::master())
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Attempted to calculate remoteWeights for a slave processor. "
            << "This is not allowed."
            << abort(FatalError);
    }
    else if (!remoteWeightsPtr_)
    {
        FatalErrorIn
        (
            "const oversetInterpolation::ListScalarFieldField&\n"
            "averageValueInterpolation::remoteWeights() const"
        )   << "Calculation of remoteWeights not possible because the data \n"
            << "exists only on the master processor. Please calculate \n"
            << "localWeights first (call .localWeights() member function)."
            << abort(FatalError);
    }

    return *remoteWeightsPtr_;
}


void Foam::leastSquaresGradientInterpolation::update() const
{
    Info<< "leastSquaresGradientInterpolation::update()" << endl;

    clearWeights();

    // Clear remote acceptor cell centres in the base class
    oversetInterpolation::update();
}


// ************************************************************************* //
