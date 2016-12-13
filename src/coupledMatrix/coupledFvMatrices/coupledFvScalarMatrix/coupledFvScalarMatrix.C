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
     Finite-Volume scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "coupledFvScalarMatrix.H"
#include "zeroGradientFvPatchFields.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
coupledSolverPerformance coupledFvMatrix<scalar>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        InfoIn("coupledFvMatrix<Type>::solve(const dictionary)")
            << "solving coupledFvMatrix<Type>" << endl;
    }

    typedef FieldField<Field, scalar> scalarFieldField;

    PtrList<lduMatrix>& matrices = *this;

    // Make a copy of the diagonal and complete the source
    scalarFieldField saveDiag(this->size());
    scalarFieldField psi(this->size());
    FieldField<Field, scalar> source(this->size());
    lduInterfaceFieldPtrsListList interfaces(this->size());

    PtrList<FieldField<Field, scalar> > bouCoeffs(this->size());
    PtrList<FieldField<Field, scalar> > intCoeffs(this->size());

    // Prepare block solution
    forAll (matrices, rowI)
    {
        fvScalarMatrix& curMatrix =
            static_cast<fvScalarMatrix&>(matrices[rowI]);

        curMatrix.completeAssembly();

        saveDiag.set(rowI, new scalarField(curMatrix.diag()));

        // Need to be able to compare references to support hacks
        // such as in jumpCyclic.  HR 17/Feb/2013
        scalarField& psiRef =
            const_cast<scalarField&>(curMatrix.psi().internalField());

        psi.set(rowI, &psiRef);
        source.set(rowI, new scalarField(curMatrix.source()));

        curMatrix.addBoundarySource(source[rowI], 0);

        interfaces[rowI] = curMatrix.psi().boundaryField().interfaces();

        curMatrix.addBoundaryDiag(curMatrix.diag(), 0);

        bouCoeffs.set
        (
            rowI,
            new FieldField<Field, scalar>
            (
                curMatrix.boundaryCoeffs()
            )
        );

        intCoeffs.set
        (
            rowI,
            new FieldField<Field, scalar>
            (
                curMatrix.internalCoeffs().component(0)
            )
        );
    }

    // Solver call
    coupledSolverPerformance solverPerf = coupledLduSolver::New
    (
        this->coupledPsiName(),
        *this,
        bouCoeffs,
        intCoeffs,
        interfaces,
        solverControls
    )->solve(psi, source, 0);

    solverPerf.print();

    // Correct boundary conditions
    forAll (matrices, rowI)
    {
        const fvScalarMatrix& curMatrix =
            static_cast<const fvScalarMatrix&>(matrices[rowI]);

        volScalarField& psiRef =
            const_cast<volScalarField&>(curMatrix.psi());

        psiRef.correctBoundaryConditions();
    }

    // Clear references to internal field without deleting the objects
    // HR 17.2.2013
    forAll (matrices, rowI)
    {
        psi.set(rowI, NULL).ptr();
    }

    return solverPerf;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
