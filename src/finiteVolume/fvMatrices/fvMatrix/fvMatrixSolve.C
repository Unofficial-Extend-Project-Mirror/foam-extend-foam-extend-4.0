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

\*---------------------------------------------------------------------------*/

#include "profiling.H"

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
Foam::lduSolverPerformance Foam::fvMatrix<Type>::solve
(
    const dictionary& solverControls
)
{
    profilingTrigger profSolve("fvMatrix::solve_" + psi_.name());

    if (debug)
    {
        Info<< "fvMatrix<Type>::solve(const dictionary&) : "
               "solving fvMatrix<Type>"
            << endl;
    }

    // Complete matrix assembly.  HJ, 17/Apr/2012
    this->completeAssembly();

    lduSolverPerformance solverPerfVec
    (
        "fvMatrix<Type>::solve",
        psi_.name()
    );

    scalarField saveDiag = diag();

    Field<Type> source = source_;

    // At this point include the boundary source from the coupled boundaries.
    // This is corrected for the implicit part by correctCmptBoundarySource
    // within the component loop.

    // Note: this is related to non-parallel coupled implicit boundaries
    // such as cyclic or cyclicGGI, which include transformation.
    // See also correctImplicitBoundarySource below.
    // HJ, 31/May/2013
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

    // Cast into a non-const to solve.  HJ, 6/May/2016
    GeometricField<Type, fvPatchField, volMesh>& psi =
       const_cast<GeometricField<Type, fvPatchField, volMesh>&>(psi_);

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

        // Correct component boundary source for the explicit part of the
        // coupled boundary conditions.  At the moment, the whole
        // coefficient-field product hass been added into the source,
        // but the implicit part for the current element needs to be taken out
        // (because it is implicit).
        // HJ, 31/May/2013
        correctImplicitBoundarySource
        (
            bouCoeffsCmpt,
            sourceCmpt,
            cmpt
        );

        lduSolverPerformance solverPerf;

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

        psi.internalField().replace(cmpt, psiCmpt);
        diag() = saveDiag;
    }

    psi.correctBoundaryConditions();

    psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerfVec);

    return solverPerfVec;
}


template<class Type>
Foam::lduSolverPerformance Foam::fvMatrix<Type>::solve()
{
    return solve(psi_.mesh().solutionDict().solverDict(psi_.name()));
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::fvMatrix<Type>::residual() const
{
    // Complete matrix assembly.  HJ, 17/Apr/2012
    fvMatrix& m = const_cast<fvMatrix&>(*this);
    m.completeAssembly();

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
