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

#include "donorSuitability.H"
#include "overlapFringe.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

defineTypeNameAndDebug(donorSuitability, 0);
defineRunTimeSelectionTable(donorSuitability, dictionary);

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::donorSuitability::donorSuitability::donorSuitability
(
    const overlapFringe& overlapFringeAlgorithm,
    const dictionary& dict
)
:
    overlapFringe_(overlapFringeAlgorithm),
    coeffDict_
    (
        dict.subDict("donorSuitability")
    ),
    dsf_(Pstream::nProcs()),
    threshold_(readScalar(coeffDict_.lookup("threshold")))
{
    // Sanity check
    if (threshold_ < SMALL)
    {
        FatalIOErrorIn
        (
            "donorSuitability::"
            "patchDistance::patchDistance()",
            coeffDict()
        )   << "Negative threshold specified. This is not allowed"
            << abort(FatalIOError);
    }
}


// * * * * * * * * * * * * * * * * Destructor* * * * * * * * * * * * * * * * //

Foam::donorSuitability::donorSuitability::~donorSuitability()
{}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

void Foam::donorSuitability::donorSuitability::combineDonorSuitabilityFunction
(
    const scalarField& localDsf
)
{
    // Perform gather-scatter for parallel run
    if (Pstream::parRun())
    {
        // Copy local donor suitability function into its processor slot
        dsf_[Pstream::myProcNo()] = localDsf;

        // Gather-scatter donor suitability function
        Pstream::gatherList(dsf_);
        Pstream::scatterList(dsf_);
    }
    // Serial run
    else
    {
        dsf_[0] = localDsf;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::overlapFringe&
Foam::donorSuitability::donorSuitability::overlapFringeAlgorithm() const
{
    return overlapFringe_;
}


const Foam::dictionary&
Foam::donorSuitability::donorSuitability::coeffDict() const
{
    return coeffDict_;
}


// ************************************************************************* //
