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
    Tetrahedral Finite Element matrix basic solvers.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Solvers * * * * * * * * * * * * * * * * * //

template<class Type>
lduSolverPerformance tetFemMatrix<Type>::solve
(
     const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve(const dictionary&) : "
               "solving tetFemMatrix<Type>"
            << endl;
    }

    // Check the matrix
    if (debug > 1)
    {
        this->check();
    }

    lduSolverPerformance solverPerfVec
    (
        "tetFemMatrix<Type>::solve",
        psi_.name()
    );

    // Add boundary source for gradient-type conditions
    addBoundarySourceDiag();

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh()().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        scalarField sourceCmpt = source_.component(cmpt);

        // Set component boundary conditions
        setComponentBoundaryConditions(cmpt, psiCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();

        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs
        (
            psi_.boundaryField().size()
        );

        FieldField<Field, scalar> coupledIntCoeffs
        (
            psi_.boundaryField().size()
        );

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

            coupledBouCoeffs.set
            (
                patchI,
                ptf.cutBouCoeffs(*this)
            );

            coupledIntCoeffs.set
            (
                patchI,
                ptf.cutIntCoeffs(*this)
            );
        }

        eliminateCouplingCoeffs();

        scalarField res(psi_.size(), 0);
        lduMatrix::residual
        (
            res,
            psiCmpt,
            sourceCmpt,
            coupledBouCoeffs,
            interfaces,
            cmpt
        );

        lduSolverPerformance solverPerf = lduSolver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            coupledBouCoeffs,
            coupledIntCoeffs,
            interfaces,
            solverControls
        )->solve(psiCmpt, sourceCmpt, cmpt);

        solverPerf.print();

        if
        (
            solverPerf.initialResidual() > solverPerfVec.initialResidual()
         && !solverPerf.singular()
        )
        {
            solverPerfVec = solverPerf;
        }

        psi_.internalField().replace(cmpt, psiCmpt);

        reconstructMatrix();
    }

    if (debug)
    {
        Info<< "tetFemMatrix<Type>::solve : correcting boundary conditions"
            << endl;
    }

    psi_.correctBoundaryConditions();

    return solverPerfVec;
}


template<class Type>
lduSolverPerformance tetFemMatrix<Type>::solve()
{
    return solve(psi_.mesh().solver(psi_.name()));
}


// Return the matrix residual
template<class Type>
tmp<Field<Type> > tetFemMatrix<Type>::residual()
{
    tmp<Field<Type> > tres(psi_.size());

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField PsiInternalCmpt = psi_.internalField().component(cmpt);
        scalarField sourceCmpt = source_.component(cmpt);

        setComponentBoundaryConditions(cmpt, PsiInternalCmpt, sourceCmpt);

        // Add the coupling coefficients
        addCouplingCoeffs();
        addCouplingSource(sourceCmpt);

        // Prepare for coupled interface update
        FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());

        forAll(psi_.boundaryField(), patchI)
        {
            const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];
            coupledBouCoeffs.set
            (
                patchI,
                new scalarField(ptf.cutBouCoeffs(*this))
            );
        }

        eliminateCouplingCoeffs();

        tres().replace
        (
            cmpt,
            lduMatrix::residual
            (
                psi_.internalField().component(cmpt),
                sourceCmpt,
                coupledBouCoeffs,
                interfaces,
                cmpt
            )
        );

        reconstructMatrix();
    }

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
