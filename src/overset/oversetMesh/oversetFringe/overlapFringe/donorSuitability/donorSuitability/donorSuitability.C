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

#include "donorSuitability.H"
#include "oversetFringe.H"

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
    const oversetFringe& oversetFringeAlgorithm,
    const dictionary& dict
)
:
    oversetFringe_(oversetFringeAlgorithm),
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


// ************************************************************************* //
