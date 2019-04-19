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

#include "loadBalanceFvMesh.H"
#include "domainDecomposition.H"
#include "fvFieldDecomposer.H"
#include "processorMeshesReconstructor.H"
#include "fvFieldReconstructor.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"

#include "passiveProcessorPolyPatch.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(loadBalanceFvMesh, 0);
    addToRunTimeSelectionTable(dynamicFvMesh, loadBalanceFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::loadBalanceFvMesh::checkLoadBalance
(
    const scalarField& weights
) const
{
    if (Pstream::parRun())
    {
        // Calculate local and global load
        scalar localLoad = sum(weights);

        scalar globalLoad = localLoad;

        reduce(globalLoad, sumOp<scalar>());

        globalLoad /= Pstream::nProcs();

        // Calculate imbalance as min of localLoad/globalLoad
        scalar imbalance = mag(1 - localLoad/globalLoad);

        reduce(imbalance, minOp<scalar>());

        Info<< "Global imbalance: " << imbalance << endl;

        if (imbalance < imbalanceTrigger_)
        {
            return true;
        }
    }

    // Serial run or low imbalance: no balancing possible
    return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::loadBalanceFvMesh::loadBalanceFvMesh(const IOobject& io)
:
    topoChangerFvMesh(io),
    dict_
    (
        IOdictionary
        (
            IOobject
            (
                "dynamicMeshDict",
                time().constant(),
                *this,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        ).subDict(typeName + "Coeffs")
    ),
    imbalanceTrigger_
    (
        readScalar(dict_.lookup("imbalanceTrigger"))
    )
{
    // Check imbalance trigger
    if (imbalanceTrigger_ < SMALL || imbalanceTrigger_ > 1)
    {
        WarningIn
        (
            "loadBalanceFvMesh::"
            "loadBalanceFvMesh(const IOobject& io)"
        )   << "Invalid imbalance trigger " << imbalanceTrigger_
            << " Should be between 0 and 1.  Resetting to 0.8"
            << endl;

        imbalanceTrigger_ = 0.8;
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::loadBalanceFvMesh::~loadBalanceFvMesh()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::loadBalanceFvMesh::update()
{
    if (Pstream::parRun())
    {
        // Decide when to balance here

        return loadBalance(dict_);
    }
    else
    {
        // Serial run: no balancing
        return false;
    }
}


// ************************************************************************* //
