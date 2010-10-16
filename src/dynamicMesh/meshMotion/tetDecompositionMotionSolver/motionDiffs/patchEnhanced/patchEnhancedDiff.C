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
    Patch-enhanced motion diffusion.

\*---------------------------------------------------------------------------*/

#include "patchEnhancedDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "elementFields.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(patchEnhancedDiff, 0);
    addToRunTimeSelectionTable(motionDiff, patchEnhancedDiff, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::patchEnhancedDiff::patchEnhancedDiff
(
    const tetDecompositionMotionSolver& mSolver
)
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
    enhance(motionGamma_);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchEnhancedDiff::~patchEnhancedDiff()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchEnhancedDiff::enhance(elementScalarField& g) const
{
    const polyMesh& m = tetMesh()();
    const polyBoundaryMesh& bdry = m.boundaryMesh();

    forAll (patchNames_, i)
    {
        label pID = bdry.findPatchID(patchNames_[i]);

        if (pID > -1)
        {
            // Cannot use patch operations: they are made for point fields
            const unallocLabelList& fc =
                m.boundaryMesh()[pID].faceCells();

            forAll (fc, fcI)
            {
                g[fc[fcI]] *= 2;
            }
        }
    }
}


// ************************************************************************* //
