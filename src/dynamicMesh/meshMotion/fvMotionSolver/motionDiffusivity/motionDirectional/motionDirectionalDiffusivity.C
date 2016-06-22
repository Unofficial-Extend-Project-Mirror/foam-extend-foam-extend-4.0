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

\*---------------------------------------------------------------------------*/

#include "motionDirectionalDiffusivity.H"
#include "surfaceInterpolate.H"
#include "zeroGradientFvPatchFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(motionDirectionalDiffusivity, 0);

    addToRunTimeSelectionTable
    (
        motionDiffusivity,
        motionDirectionalDiffusivity,
        Istream
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::motionDirectionalDiffusivity
(
    const fvMotionSolver& mSolver,
    Istream& mdData
)
:
    uniformDiffusivity(mSolver, mdData),
    diffusivityVector_(mdData)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::motionDirectionalDiffusivity::~motionDirectionalDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::motionDirectionalDiffusivity::correct()
{
    const fvMesh& mesh = mSolver().mesh();

    static bool first = true;

    if (!first)
    {
        const volVectorField& cellMotionU =
            mesh.lookupObject<volVectorField>("cellMotionU");

        volVectorField D
        (
            IOobject
            (
                "D",
                mesh.time().timeName(),
                mesh
            ),
            diffusivityVector_.y()*vector::one
          + (diffusivityVector_.x() - diffusivityVector_.y())*cellMotionU
           /(mag(cellMotionU) + dimensionedScalar("small", dimVelocity, SMALL)),
            zeroGradientFvPatchVectorField::typeName
        );
        D.correctBoundaryConditions();

        surfaceVectorField n = mesh.Sf()/mesh.magSf();
        faceDiffusivity_ == (n & cmptMultiply(fvc::interpolate(D), n));
    }
    else
    {
        first = false;
        const_cast<fvMotionSolver&>(mSolver()).solve();
        correct();
    }
}


// ************************************************************************* //
