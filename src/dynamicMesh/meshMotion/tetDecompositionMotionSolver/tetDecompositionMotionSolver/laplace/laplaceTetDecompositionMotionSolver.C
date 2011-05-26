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
    Mesh motion solver for a polyMesh.  Based on solving the
    vertex-based laplace motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "laplaceTetDecompositionMotionSolver.H"
#include "motionDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "tetFem.H"
#include "elementFields.H"
#include "tetFec.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(laplaceTetDecompositionMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        tetDecompositionMotionSolver,
        laplaceTetDecompositionMotionSolver,
        dictionary
    );
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::laplaceTetDecompositionMotionSolver::laplaceTetDecompositionMotionSolver
(
    const polyMesh& mesh,
    Istream&
)
:
    tetDecompositionMotionSolver(mesh),
    diffusionPtr_(motionDiff::New(*this).ptr()),
    firstMotion_(true),
    solverPerf_()
{
    frozen_ = Switch(lookup("frozenDiffusion"));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::laplaceTetDecompositionMotionSolver::
~laplaceTetDecompositionMotionSolver()
{
    deleteDemandDrivenData(diffusionPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::laplaceTetDecompositionMotionSolver::solve()
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


void Foam::laplaceTetDecompositionMotionSolver::updateMesh
(
    const mapPolyMesh& mpm
)
{
    firstMotion_ = true;
    tetDecompositionMotionSolver::updateMesh(mpm);
}


// ************************************************************************* //
