/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

Description
    Linear distance-based motion diffusion.

\*---------------------------------------------------------------------------*/

#include "linearDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"
#include "patchWave.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(linearDiff, 0);
    addToRunTimeSelectionTable(motionDiff, linearDiff, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::linearDiff::linearDiff(const tetMotionSolver& mSolver)
:
    motionDiff(mSolver),
    patchNames_(mSolver.lookup("distancePatches")),
    motionGamma_
    (
        IOobject
        (
            "motionGamma",
            tetMesh().time().timeName(),
            tetMesh()(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        tetMesh(),
        dimensionedScalar("1.0", dimless, 1.0)
    )
{
    motionGamma_.internalField() = 1.0/L();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::linearDiff::~linearDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField> Foam::linearDiff::L() const
{
    const polyMesh& m = tetMesh()();
    const polyBoundaryMesh& bdry = m.boundaryMesh();

    labelHashSet patchSet(bdry.size());

    forAll (patchNames_, i)
    {
        label pID = bdry.findPatchID(patchNames_[i]);

        if (pID > -1)
        {
            patchSet.insert(pID);
        }
    }

    if (patchSet.size() > 0)
    {
        return tmp<scalarField>
        (
            new scalarField(patchWave(m, patchSet, false).distance())
        );
    }
    else
    {
        return tmp<scalarField>(new scalarField(m.nCells(), 1.0));
    }
}


void Foam::linearDiff::correct()
{
    motionGamma().internalField() = 1.0/L();
}


// ************************************************************************* //
