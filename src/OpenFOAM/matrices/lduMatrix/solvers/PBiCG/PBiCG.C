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

#include "PBiCG.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(PBiCG, 0);

    lduSolver::addasymMatrixConstructorToTable<PBiCG>
        addPBiCGAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PBiCG::PBiCG
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
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduSolverPerformance Foam::PBiCG::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // --- Setup class containing solver performance data
    lduSolverPerformance solverPerf
    (
        lduMatrix::preconditioner::getName(dict()) + typeName,
        fieldName()
    );


    register label nCells = x.size();

    scalar* __restrict__ xPtr = x.begin();

    scalarField pA(nCells);
    scalar* __restrict__ pAPtr = pA.begin();

    scalarField pT(nCells, 0.0);
    scalar* __restrict__ pTPtr = pT.begin();

    scalarField wA(nCells);
    scalar* __restrict__ wAPtr = wA.begin();

    scalarField wT(nCells);
    scalar* __restrict__ wTPtr = wT.begin();

    scalar wArT = matrix_.great_;
    scalar wArTold = wArT;

    // Calculate A.x and T.x
    matrix_.Amul(wA, x, coupleBouCoeffs_, interfaces_, cmpt);
    matrix_.Tmul(wT, x, coupleIntCoeffs_, interfaces_, cmpt);

    // Calculate initial residual and transpose residual fields
    scalarField rA(b - wA);
    scalarField rT(b - wT);
    scalar* __restrict__ rAPtr = rA.begin();
    scalar* __restrict__ rTPtr = rT.begin();

    // Calculate normalisation factor
    scalar normFactor = this->normFactor(x, b, wA, pA, cmpt);

    if (lduMatrix::debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate normalised residual norm
    solverPerf.initialResidual() = gSumMag(rA)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // Check convergence, solve if not converged
    if (!stop(solverPerf))
    {
        // Select and construct the preconditioner
        autoPtr<lduPreconditioner> preconPtr;

        preconPtr =
            lduPreconditioner::New
            (
                matrix_,
                coupleBouCoeffs_,
                coupleIntCoeffs_,
                interfaces_,
                dict()
            );

        // Solver iteration
        do
        {
            // Store previous wArT
            wArTold = wArT;

            // Precondition residuals
            preconPtr->precondition(wA, rA, cmpt);
            preconPtr->preconditionT(wT, rT, cmpt);

            // Update search directions:
            wArT = gSumProd(wA, rT);

            if (solverPerf.nIterations() == 0)
            {
                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell];
                    pTPtr[cell] = wTPtr[cell];
                }
            }
            else
            {
                scalar beta = wArT/wArTold;

                for (register label cell=0; cell<nCells; cell++)
                {
                    pAPtr[cell] = wAPtr[cell] + beta*pAPtr[cell];
                    pTPtr[cell] = wTPtr[cell] + beta*pTPtr[cell];
                }
            }


            // Update preconditioned residuals
            matrix_.Amul(wA, pA, coupleBouCoeffs_, interfaces_, cmpt);
            matrix_.Tmul(wT, pT, coupleIntCoeffs_, interfaces_, cmpt);

            scalar wApT = gSumProd(wA, pT);


            // Test for singularity
            if (solverPerf.checkSingularity(mag(wApT)/normFactor)) break;


            // Update solution and residual:

            scalar alpha = wArT/wApT;

            for (register label cell=0; cell<nCells; cell++)
            {
                xPtr[cell] += alpha*pAPtr[cell];
                rAPtr[cell] -= alpha*wAPtr[cell];
                rTPtr[cell] -= alpha*wTPtr[cell];
            }

            solverPerf.finalResidual() = gSumMag(rA)/normFactor;
            solverPerf.nIterations()++;
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


// ************************************************************************* //
