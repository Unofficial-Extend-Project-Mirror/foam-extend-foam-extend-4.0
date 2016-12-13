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
    vertex-based pseudoSolid motion equation.  The boundary motion is set as a
    boundary condition on the motion velocity variable motionU.

\*---------------------------------------------------------------------------*/

#include "pseudoSolidTetMotionSolver.H"
#include "motionDiff.H"
#include "addToRunTimeSelectionTable.H"
#include "tetFem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pseudoSolidTetMotionSolver, 0);
    addToRunTimeSelectionTable
    (
        tetMotionSolver,
        pseudoSolidTetMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::pseudoSolidTetMotionSolver::
pseudoSolidTetMotionSolver
(
    const polyMesh& mesh,
    Istream& msData
)
:
    laplaceTetMotionSolver(mesh, msData)
{
    const dictionary& pseudoSolidDic = subDict("pseudoSolid");

    nu_ = readScalar(pseudoSolidDic.lookup("poissonsRatio"));

    nCorrectors_ =  readInt(pseudoSolidDic.lookup("nCorrectors"));

    convergenceTolerance_ =
        readScalar(pseudoSolidDic.lookup("convergenceTolerance"));
}



// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::pseudoSolidTetMotionSolver::
~pseudoSolidTetMotionSolver()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pseudoSolidTetMotionSolver::solve()
{
    // Solve for mesh motion

    if (!frozen_ && !firstMotion_)
    {
        Pout << "Correct mesh motion diffusion field." << endl;

        diffusionPtr_->correct();
    }

    int iCorr = 0;
    scalar initialResidual = 0;

    do
    {
        Pout << "Correction: " << ++iCorr << endl;

        tetPointVectorField& U = motionU();
        const elementScalarField& mu = diffusion().motionGamma();

        tetFemVectorMatrix motionEqn
        (
            tetFem::laplacian(mu, U)
          + tetFem::laplacianTranspose(mu, U)
          + tetFem::laplacianTrace
            (
                (2*nu_/(1 - 2*nu_))*mu,
                motionU()
            )
        );

        // Apply motion constraints
        applyConstraints(motionEqn);

        // Solve the motion equation
        initialResidual = motionEqn.solve().initialResidual();

        Pout << "Initial residual: " << initialResidual << endl;
    }
    while (initialResidual > convergenceTolerance_ && iCorr < nCorrectors_);
}


// ************************************************************************* //
