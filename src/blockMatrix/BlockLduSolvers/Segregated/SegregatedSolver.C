/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Description
    Segregated solver for symmetric and asymmetric matrices.
    Calls scalar solver component-by-component

\*---------------------------------------------------------------------------*/

#include "SegregatedSolver.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from matrix and solver data stream
template<class Type>
Foam::SegregatedSolver<Type>::SegregatedSolver
(
    const word& fieldName,
    const BlockLduMatrix<Type>& matrix,
    const dictionary& dict
)
:
    BlockLduSolver<Type>(fieldName, matrix, dict),
    scalarX_(matrix.lduAddr().size()),
    scalarMatrix_(matrix.mesh()),
    scalarB_(matrix.lduAddr().size())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::BlockSolverPerformance<Type> Foam::SegregatedSolver<Type>::solve
(
    Field<Type>& x,
    const Field<Type>& b
)
{
    // Get reference to matrix, x and b
    const BlockLduMatrix<Type>& blockMatrix = this->matrix_;

    // Determine if the diagonal or off-diagonal is expanded

    // Check switching diagonal
    bool switchingDiag = false;

    if (blockMatrix.diag().activeType() == blockCoeffBase::SCALAR)
    {
        scalarMatrix_.diag() = blockMatrix.diag().asScalar();
    }
    else
    {
        switchingDiag = true;
    }

    // Check switching upper
    bool switchingUpper = false;

    if (blockMatrix.thereIsUpper())
    {
        if (blockMatrix.upper().activeType() == blockCoeffBase::SCALAR)
        {
            scalarMatrix_.upper() = blockMatrix.upper().asScalar();
        }
        else
        {
            switchingUpper = true;
        }
    }

    // Check switching lower
    bool switchingLower = false;

    if (blockMatrix.thereIsLower())
    {
        if (blockMatrix.lower().activeType() == blockCoeffBase::SCALAR)
        {
            scalarMatrix_.lower() = blockMatrix.lower().asScalar();
        }
        else
        {
            switchingLower = true;
        }
    }

    // Decouple quadratic coupling by multiplying out the square coefficient
    // coupling
    Field<Type>* dBPtr = NULL;

    if (blockMatrix.componentCoupled())
    {
        if (BlockLduMatrix<Type>::debug >= 2)
        {
            Info << " Component coupled segregation" << endl;
        }

        dBPtr = new Field<Type>(b);
        blockMatrix.segregateB(*dBPtr, x);
    }

    // Prepare solver performance
    word segSolverName(this->dict().lookup("solver"));

    BlockSolverPerformance<Type> solverPerf
    (
        typeName + "_" + segSolverName,
        this->fieldName()
    );

    // Loop through the form component by component and call the
    // scalar solver
    for
    (
        direction cmpt = 0;
        cmpt < pTraits<Type>::nComponents;
        cmpt++
    )
    {
        // Grab the x and b for the current component
        scalarX_ = x.component(cmpt);

        // Handle b decoupling
        if (dBPtr)
        {
            scalarB_ = dBPtr->component(cmpt);
        }
        else
        {
            scalarB_ = b.component(cmpt);
        }

        // Switch diagonal, upper and lower
        if (switchingDiag)
        {
            scalarMatrix_.diag() = blockMatrix.diag().component(cmpt);
        }

        if (switchingUpper)
        {
            scalarMatrix_.upper() = blockMatrix.upper().component(cmpt);
        }

        if (switchingLower)
        {
            scalarMatrix_.lower() = blockMatrix.lower().component(cmpt);
        }

        // Call the scalar solver
        BlockSolverPerformance<scalar> scalarPerf =
            blockScalarSolver::New
            (
                this->fieldName(),
                scalarMatrix_,
                this->dict()
            )->solve(scalarX_, scalarB_);

        // Replace the solution component
        x.replace(cmpt, scalarX_);

        // Grab solver performance and combine
        solverPerf.initialResidual().replace
        (
            cmpt,
            scalarPerf.initialResidual()
        );

        solverPerf.finalResidual().replace
        (
            cmpt,
            scalarPerf.finalResidual()
        );

        solverPerf.nIterations() =
            max
            (
                solverPerf.nIterations(),
                scalarPerf.nIterations()
            );

        solverPerf.converged() =
            (solverPerf.converged() && scalarPerf.converged());

        solverPerf.singular() =
            solverPerf.singular() && scalarPerf.singular();
    }

    // Clear decoupled b if allocated
    if (dBPtr)
    {
        delete dBPtr;
        dBPtr = NULL;
    }

    return solverPerf;
}


// ************************************************************************* //
