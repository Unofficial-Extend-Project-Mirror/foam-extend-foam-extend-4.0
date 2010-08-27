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

\*---------------------------------------------------------------------------*/

#include "smoothSolver.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(smoothSolver, 0);

    lduSolver::addsymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverSymMatrixConstructorToTable_;

    lduSolver::addasymMatrixConstructorToTable<smoothSolver>
        addsmoothSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::smoothSolver::smoothSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& dict
)
:
    lduSolver
    (
        fieldName,
        matrix,
        coupleBouCoeffs,
        coupleIntCoeffs,
        interfaces,
        dict
    )
{
    readControls();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::smoothSolver::readControls()
{
    lduSolver::readControls();
    nSweeps_ = dict().lookupOrDefault<label>("nSweeps", 1);
}


Foam::lduSolverPerformance Foam::smoothSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    lduSolverPerformance solverPerf(typeName, fieldName());

    // If the nSweeps_ is negative do a fixed number of sweeps
    if (nSweeps_ < 0)
    {
        autoPtr<lduMatrix::smoother> smootherPtr = lduMatrix::smoother::New
        (
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaces_,
            dict()
        );

        smootherPtr->smooth
        (
            x,
            b,
            cmpt,
            minIter()
        );

        solverPerf.nIterations() += minIter();
    }

    // HJ, bug fix.  Now do normal sweeps.  HJ, 19/Jan/2009
    scalar normFactor = 0;

    {
        scalarField Ax(x.size());
        scalarField temp(x.size());

        // Calculate A.x
        matrix_.Amul(Ax, x, coupleBouCoeffs_, interfaces_, cmpt);

        // Calculate normalisation factor
        normFactor = this->normFactor(x, b, Ax, temp, cmpt);

        // Calculate residual magnitude
        solverPerf.initialResidual() = gSumMag(b - Ax)/normFactor;
        solverPerf.finalResidual() = solverPerf.initialResidual();
    }

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }


    // Check convergence, solve if not converged
    if (!stop(solverPerf))
    {
        autoPtr<lduMatrix::smoother> smootherPtr =
            lduMatrix::smoother::New
            (
                matrix_,
                coupleBouCoeffs_,
                coupleIntCoeffs_,
                interfaces_,
                dict()
            );

        // Smoothing loop
        do
        {
            smootherPtr->smooth
            (
                x,
                b,
                cmpt,
                nSweeps_
            );

            // Calculate the residual to check convergence
            solverPerf.finalResidual() = gSumMag
            (
                matrix_.residual
                (
                    x,
                    b,
                    coupleBouCoeffs_,
                    interfaces_,
                    cmpt
                )
            )/normFactor;

            solverPerf.nIterations() += nSweeps_;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
