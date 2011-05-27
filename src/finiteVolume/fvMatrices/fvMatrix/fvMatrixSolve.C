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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvMatrix<Type>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction cmpt,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei].component(cmpt) +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
               *value;
        }
    }
}


template<class Type>
Foam::lduMatrix::solverPerformance Foam::fvMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "fvMatrix<Type>::solve(const dictionary&) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    lduSolverPerformance solverPerfVec
    (
        "fvMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implicit part by updateMatrixInterfaces within
    // the component loop.
    addBoundarySource(source);

    typename Type::labelType validComponents
    (
        pow
        (
            psi_.mesh().solutionD(),
            pTraits<typename powProduct<Vector<label>, Type::rank>::type>::zero
        )
    );

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        if (validComponents[cmpt] == -1) continue;

        // Copy field and source

        scalarField psiCmpt = psi_.internalField().component(cmpt);
        addBoundaryDiag(diag(), cmpt);

        scalarField sourceCmpt = source.component(cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        FieldField<Field, scalar> intCoeffsCmpt
        (
            internalCoeffs_.component(cmpt)
        );

        // Use the initMatrixInterfaces and updateMatrixInterfaces to correct
        // bouCoeffsCmpt for the explicit part of the coupled boundary
        // conditions
        initMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        lduMatrix::solverPerformance solverPerf;

        // Solver call
        solverPerf = lduMatrix::solver::New
        (
            psi_.name() + pTraits<Type>::componentNames[cmpt],
            *this,
            bouCoeffsCmpt,
            intCoeffsCmpt,
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
        diag() = saveDiag;
    }

    psi_.correctBoundaryConditions();

    return solverPerfVec;
}


template<class Type>
Foam::autoPtr<typename Foam::fvMatrix<Type>::fvSolver>
Foam::fvMatrix<Type>::solver()
{
    return solver(psi_.mesh().solverDict(psi_.name()));
}

template<class Type>
Foam::lduMatrix::solverPerformance Foam::fvMatrix<Type>::fvSolver::solve()
{
    return solve(psi_.mesh().solverDict(psi_.name()));
}


template<class Type>
Foam::lduMatrix::solverPerformance Foam::fvMatrix<Type>::solve()
{
    return solve(psi_.mesh().solverDict(psi_.name()));
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fvMatrix<Type>::residual() const
{
    // Bug fix: Creating a tmp out of a const reference will change the field
    // HJ, 15/Apr/2011
    tmp<Field<Type> > tres(new Field<Type>(source_));
    Field<Type>& res = tres();

    addBoundarySource(res);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Loop over field components
    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        scalarField psiCmpt = psi_.internalField().component(cmpt);

        scalarField boundaryDiagCmpt(psi_.size(), 0);
        addBoundaryDiag(boundaryDiagCmpt, cmpt);

        FieldField<Field, scalar> bouCoeffsCmpt
        (
            boundaryCoeffs_.component(cmpt)
        );

        res.replace
        (
            cmpt,
            lduMatrix::residual
            (
                psiCmpt,
                res.component(cmpt) - boundaryDiagCmpt*psiCmpt,
                bouCoeffsCmpt,
                interfaces,
                cmpt
            )
        );
    }

    return tres;
}


// ************************************************************************* //
