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

#include "GAMGSolver.H"
#include "ICCG.H"
#include "BICCG.H"
#include "SubField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::lduMatrix::solverPerformance Foam::GAMGSolver::solve
(
    scalarField& x,
    const scalarField& b,
    const direction cmpt
) const
{
    // Setup class containing solver performance data
    lduSolverPerformance solverPerf(typeName, fieldName());

    // Calculate A.x used to calculate the initial residual
    scalarField Ax(x.size());
    matrix_.Amul(Ax, x, coupleBouCoeffs_, interfaces_, cmpt);

    // Create the storage for the finestCorrection which may be used as a
    // temporary in normFactor
    scalarField finestCorrection(x.size());

    // Calculate normalisation factor
    scalar normFactor = this->normFactor(x, b, Ax, finestCorrection, cmpt);

    if (debug >= 2)
    {
        Info<< "   Normalisation factor = " << normFactor << endl;
    }

    // Calculate initial finest-grid residual field
    scalarField finestResidual(b - Ax);

    // Calculate normalised residual for convergence test
    solverPerf.initialResidual() = gSumMag(finestResidual)/normFactor;
    solverPerf.finalResidual() = solverPerf.initialResidual();


    // Check convergence, solve if not converged
    if (!stop(solverPerf))
    {
        // Create coarse grid correction fields
        PtrList<scalarField> coarseCorrX;

        // Create coarse grid bs
        PtrList<scalarField> coarseB;

        // Create the smoothers for all levels
        PtrList<lduSmoother> smoothers;

        // Initialise the above data structures
        initVcycle(coarseCorrX, coarseB, smoothers);

        do
        {
            Vcycle
            (
                smoothers,
                x,
                b,
                Ax,
                finestCorrection,
                finestResidual,
                coarseCorrX,
                coarseB,
                cmpt
            );

            // Calculate finest level residual field
            matrix_.Amul(Ax, x, coupleBouCoeffs_, interfaces_, cmpt);
            finestResidual = b;
            finestResidual -= Ax;

            solverPerf.finalResidual() = gSumMag(finestResidual)/normFactor;
            solverPerf.nIterations()++;
            if (debug >= 2)
            {
                solverPerf.print();
            }
        } while (!stop(solverPerf));
    }

    return solverPerf;
}


void Foam::GAMGSolver::Vcycle
(
    const PtrList<lduSmoother>& smoothers,
    scalarField& x,
    const scalarField& b,
    scalarField& Ax,
    scalarField& finestCorrection,
    scalarField& finestResidual,
    PtrList<scalarField>& coarseCorrX,
    PtrList<scalarField>& coarseB,
    const direction cmpt
) const
{
    //debug = 2;

    const label coarsestLevel = matrixLevels_.size() - 1;

    // Restrict finest grid residual for the next level up
    agglomeration_.restrictField(coarseB[0], finestResidual, 0);

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< "Pre-smoothing scaling factors: ";
    }


    // Residual restriction (going to coarser levels)
    for (label leveli = 0; leveli < coarsestLevel; leveli++)
    {
        // If the optional pre-smoothing sweeps are selected
        // smooth the coarse-grid field for the restricted b
        if (nPreSweeps_)
        {
            coarseCorrX[leveli] = 0.0;

            smoothers[leveli + 1].smooth
            (
                coarseCorrX[leveli],
                coarseB[leveli],
                cmpt,
                nPreSweeps_ + leveli
            );

            scalarField::subField ACf
            (
                Ax,
                coarseCorrX[leveli].size()
            );

            // Scale coarse-grid correction field
            // but not on the coarsest level because it evaluates to 1
            if (scaleCorrection_ && leveli < coarsestLevel - 1)
            {
                scalar sf = scalingFactor
                (
                    const_cast<scalarField&>(ACf.operator const scalarField&()),
                    matrixLevels_[leveli],
                    coarseCorrX[leveli],
                    coupleLevelsBouCoeffs_[leveli],
                    interfaceLevels_[leveli],
                    coarseB[leveli],
                    cmpt
                );

                if (debug >= 2)
                {
                    Pout<< sf << " ";
                }

                coarseCorrX[leveli] *= sf;
            }

            // Correct the residual with the new solution
            matrixLevels_[leveli].Amul
            (
                const_cast<scalarField&>(ACf.operator const scalarField&()),
                coarseCorrX[leveli],
                coupleLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                cmpt
            );

            coarseB[leveli] -= ACf;
        }

        // Residual is equal to b
        agglomeration_.restrictField
        (
            coarseB[leveli + 1],
            coarseB[leveli],
            leveli + 1
        );
    }

    if (debug >= 2 && nPreSweeps_)
    {
        Pout<< endl;
    }


    // Solve Coarsest level with either an iterative or direct solver
    solveCoarsestLevel
    (
        coarseCorrX[coarsestLevel],
        coarseB[coarsestLevel]
    );


    if (debug >= 2)
    {
        Pout<< "Post-smoothing scaling factors: ";
    }

    // Smoothing and prolongation of the coarse correction fields
    // (going to finer levels)
    for (label leveli = coarsestLevel - 1; leveli >= 0; leveli--)
    {
        // Create a field for the pre-smoothed correction field
        // as a sub-field of the finestCorrection which is not
        // currently being used
        scalarField::subField preSmoothedCoarseCorrField
        (
            finestCorrection,
            coarseCorrX[leveli].size()
        );

        // Only store the preSmoothedCoarseCorrField is pre-smoothing is used
        if (nPreSweeps_)
        {
            preSmoothedCoarseCorrField.assign(coarseCorrX[leveli]);
        }

        agglomeration_.prolongField
        (
            coarseCorrX[leveli],
            coarseCorrX[leveli + 1],
            leveli + 1
        );

        // Scale coarse-grid correction field
        // but not on the coarsest level because it evaluates to 1
        if (scaleCorrection_ && leveli < coarsestLevel - 1)
        {
            // Create A.x for this coarse level as a sub-field of Ax
            scalarField::subField ACf
            (
                Ax,
                coarseCorrX[leveli].size()
            );

            scalar sf = scalingFactor
            (
                const_cast<scalarField&>(ACf.operator const scalarField&()),
                matrixLevels_[leveli],
                coarseCorrX[leveli],
                coupleLevelsBouCoeffs_[leveli],
                interfaceLevels_[leveli],
                coarseB[leveli],
                cmpt
            );


            if (debug >= 2)
            {
                Pout<< sf << " ";
            }

            coarseCorrX[leveli] *= sf;
        }

        // Only add the preSmoothedCoarseCorrField is pre-smoothing is used
        if (nPreSweeps_)
        {
            coarseCorrX[leveli] += preSmoothedCoarseCorrField;
        }

        smoothers[leveli + 1].smooth
        (
            coarseCorrX[leveli],
            coarseB[leveli],
            cmpt,
            nPostSweeps_ + leveli
        );
    }

    // Prolong the finest level correction
    agglomeration_.prolongField
    (
        finestCorrection,
        coarseCorrX[0],
        0
    );

    if (scaleCorrection_)
    {
        // Calculate finest level scaling factor
        scalar fsf = scalingFactor
        (
            Ax,
            matrix_,
            finestCorrection,
            coupleBouCoeffs_,
            interfaces_,
            finestResidual,
            cmpt
        );

        if (debug >= 2)
        {
            Pout<< fsf << endl;
        }

        forAll(x, i)
        {
            x[i] += fsf*finestCorrection[i];
        }
    }
    else
    {
        forAll(x, i)
        {
            x[i] += finestCorrection[i];
        }
    }

    smoothers[0].smooth
    (
        x,
        b,
        cmpt,
        nFinestSweeps_
    );
}


void Foam::GAMGSolver::initVcycle
(
    PtrList<scalarField>& coarseCorrX,
    PtrList<scalarField>& coarseB,
    PtrList<lduSmoother>& smoothers
) const
{
    coarseCorrX.setSize(matrixLevels_.size());
    coarseB.setSize(matrixLevels_.size());
    smoothers.setSize(matrixLevels_.size() + 1);

    // Create the smoother for the finest level
    smoothers.set
    (
        0,
        lduSmoother::New
        (
            matrix_,
            coupleBouCoeffs_,
            coupleIntCoeffs_,
            interfaces_,
            dict()
        )
    );

    forAll (matrixLevels_, leveli)
    {
        coarseCorrX.set
        (
            leveli,
            new scalarField
            (
                agglomeration_.meshLevel(leveli + 1).lduAddr().size()
            )
        );

        coarseB.set
        (
            leveli,
            new scalarField
            (
                agglomeration_.meshLevel(leveli + 1).lduAddr().size()
            )
        );

        smoothers.set
        (
            leveli + 1,
            lduSmoother::New
            (
                matrixLevels_[leveli],
                coupleLevelsBouCoeffs_[leveli],
                coupleLevelsIntCoeffs_[leveli],
                interfaceLevels_[leveli],
                dict()
            )
        );
    }
}


void Foam::GAMGSolver::solveCoarsestLevel
(
    scalarField& coarsestCorrX,
    const scalarField& coarsestB
) const
{
    if (directSolveCoarsest_)
    {
        coarsestCorrX = coarsestB;
        coarsestLUMatrixPtr_->solve(coarsestCorrX);
    }
    else
    {
        const label coarsestLevel = matrixLevels_.size() - 1;
        coarsestCorrX = 0;
        lduSolverPerformance coarseSolverPerf;

        if (matrixLevels_[coarsestLevel].asymmetric())
        {
            coarseSolverPerf = BICCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                coupleLevelsBouCoeffs_[coarsestLevel],
                coupleLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance(),
                relTolerance()
            ).solve
            (
                coarsestCorrX,
                coarsestB
            );
        }
        else
        {
            coarseSolverPerf = ICCG
            (
                "coarsestLevelCorr",
                matrixLevels_[coarsestLevel],
                coupleLevelsBouCoeffs_[coarsestLevel],
                coupleLevelsIntCoeffs_[coarsestLevel],
                interfaceLevels_[coarsestLevel],
                tolerance(),
                relTolerance()
            ).solve
            (
                coarsestCorrX,
                coarsestB
            );
        }

        if (debug >= 2)
        {
            coarseSolverPerf.print();
        }
    }
}


// ************************************************************************* //
