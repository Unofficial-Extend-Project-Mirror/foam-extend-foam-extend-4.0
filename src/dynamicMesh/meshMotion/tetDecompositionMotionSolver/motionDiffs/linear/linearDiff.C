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
Foam::linearDiff::linearDiff(const tetDecompositionMotionSolver& mSolver)
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
        return
            tmp<scalarField>
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
