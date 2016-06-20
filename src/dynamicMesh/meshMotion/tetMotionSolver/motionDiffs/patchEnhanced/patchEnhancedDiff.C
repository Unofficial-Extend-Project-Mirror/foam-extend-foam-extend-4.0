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
    const tetMotionSolver& mSolver
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
