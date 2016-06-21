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
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based laplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "laplaceTetMotionSolver.H"
#include "motionDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "tetFem.H"
#include "elementFields.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplaceTetMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        tetMotionSolver,
        laplaceTetMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::laplaceTetMotionSolver::laplaceTetMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    tetMotionSolver(mesh),
    diffusionPtr_(motionDiff::New(*this).ptr()),
    firstMotion_(true),
    solverPerf_()
{
    frozen_ = Switch(lookup("frozenDiffusion"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laplaceTetMotionSolver::
~laplaceTetMotionSolver()
{
    deleteDemandDrivenData(diffusionPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laplaceTetMotionSolver::solve()
{
    // Solve for mesh motion

    if (!frozen_ && !firstMotion_)
    {
        Info << "Correct mesh motion diffusion field." << endl;

        diffusionPtr_->correct();
    }

    tetFemVectorMatrix motionEqn
    (
        tetFem::laplacian
        (
            diffusion().motionGamma(),
            motionU()
        )
    );

    // Apply motion constraints
    applyConstraints(motionEqn);

    // Solve the motion equation
    if (firstMotion_)
    {
        firstMotion_ = false;

        // In the first solution, solve the motion twice to avoid relative
        // tolerance problem
        for (label i = 0; i < 2; i++)
        {
            solverPerf_ = motionEqn.solve();
        }
    }
    else
    {
        solverPerf_ = motionEqn.solve();
    }


    if (needTotDisplacement())
    {
        totDisplacement() += motionU()*tetMesh().time().deltaT();
    }
}


void Foam::laplaceTetMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    firstMotion_ = true;
    tetMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
