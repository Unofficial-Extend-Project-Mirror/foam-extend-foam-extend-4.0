/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-7 H. Jasak All rights reserved
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
    Coupled Finite Volume matrix

\*---------------------------------------------------------------------------*/

#include "coupledLduSolver.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::coupledFvMatrix<Type>::checkSize() const
{
    if (this->size() < 1)
    {
        FatalErrorIn("void Foam::coupledFvMatrix<Type>::checkSize() const")
            << "No matrices added to coupled matrix"
            << abort(FatalError);
    }
}


template<class Type>
Foam::word Foam::coupledFvMatrix<Type>::coupledPsiName() const
{
    word cpn;

    const PtrList<lduMatrix>& matrices = *this;

    forAll (matrices, rowI)
    {
        const fvMatrix<Type>& curMatrix =
            static_cast<const fvMatrix<Type>& >(matrices[rowI]);

        cpn += curMatrix.psi().name();

        if (rowI < matrices.size() - 1)
        {
            cpn += "+";
        }
    }

    return cpn;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::coupledSolverPerformance
Foam::coupledFvMatrix<Type>::solve(const dictionary& solverControls)
{
    if (debug)
    {
        InfoIn("coupledFvMatrix<Type>::solve(const dictionary&)")
            << "solving coupledFvMatrix<Type>" << endl;
    }

    coupledSolverPerformance solverPerfVec
    (
        "coupledFvMatrix<Type>::solve",
        this->coupledPsiName()
    );

    typedef FieldField<Field, scalar> scalarFieldField;

    PtrList<lduMatrix>& matrices = *this;

    // Make a copy of the diagonal and complete the source
    scalarFieldField saveDiag(this->size());
    scalarFieldField psiCmpt(this->size());
    FieldField<Field, Type> source(this->size());
    scalarFieldField sourceCmpt(this->size());
    lduInterfaceFieldPtrsListList interfaces(this->size());

    // Prepare block solution
    forAll (matrices, rowI)
    {
        fvMatrix<Type>& curMatrix =
            static_cast<fvMatrix<Type>& >(matrices[rowI]);

        saveDiag.set(rowI, new scalarField(curMatrix.diag()));
        psiCmpt.set(rowI, new scalarField(curMatrix.psi().size()));
        source.set(rowI, new Field<Type>(curMatrix.source()));
        sourceCmpt.set(rowI, new scalarField(curMatrix.psi().size()));

        curMatrix.addBoundarySource(source[rowI]);

        interfaces[rowI] = curMatrix.psi().boundaryField().interfaces();
    }

    for (direction cmpt = 0; cmpt < Type::nComponents; cmpt++)
    {
        // Copy field and source

        PtrList<FieldField<Field, scalar> > bouCoeffsCmpt(this->size());
        PtrList<FieldField<Field, scalar> > intCoeffsCmpt(this->size());

        forAll (matrices, rowI)
        {

            fvMatrix<Type>& curMatrix =
                static_cast<fvMatrix<Type>& >(matrices[rowI]);

            psiCmpt[rowI] = curMatrix.psi().internalField().component(cmpt);
            curMatrix.addBoundaryDiag(curMatrix.diag(), cmpt);

            sourceCmpt[rowI] = source[rowI].component(cmpt);

            bouCoeffsCmpt.set
            (
                rowI,
                new FieldField<Field, scalar>
                (
                    curMatrix.boundaryCoeffs().component(cmpt)
                )
            );

            intCoeffsCmpt.set
            (
                rowI,
                new FieldField<Field, scalar>
                (
                    curMatrix.internalCoeffs().component(cmpt)
                )
            );
        }

        this->initMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        this->updateMatrixInterfaces
        (
            bouCoeffsCmpt,
            interfaces,
            psiCmpt,
            sourceCmpt,
            cmpt
        );

        coupledSolverPerformance solverPerf
        (
            "coupledFvMatrix<Type>::solve",
            this->coupledPsiName()
        );

        // Solver call
        solverPerf = coupledLduSolver::New
        (
            this->coupledPsiName() + pTraits<Type>::componentNames[cmpt],
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

        // Update solution
        forAll (matrices, rowI)
        {
            fvMatrix<Type>& curMatrix =
                static_cast<fvMatrix<Type>& >(matrices[rowI]);
            curMatrix.psi().internalField().replace(cmpt, psiCmpt[rowI]);
            curMatrix.diag() = saveDiag[rowI];
        }
    }

    // Correct boundart conditions
    forAll (matrices, rowI)
    {
        fvMatrix<Type>& curMatrix =
            static_cast<fvMatrix<Type>& >(matrices[rowI]);
        curMatrix.psi().correctBoundaryConditions();
    }

    return solverPerfVec;
}


template<class Type>
Foam::coupledSolverPerformance Foam::coupledFvMatrix<Type>::solve()
{
    // Check matrices are present.  Using first matrix for controls
    checkSize();

    const fvMatrix<Type>& m =
        static_cast<const fvMatrix<Type>& >(this->operator[](0));

    return solve(m.psi().mesh().solver(coupledPsiName()));
}


// ************************************************************************* //
