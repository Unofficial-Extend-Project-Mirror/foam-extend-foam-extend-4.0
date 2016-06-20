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

#include "fvScalarMatrix.H"
#include "zeroGradientFvPatchFields.H"

#include "profiling.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
void Foam::fvMatrix<Foam::scalar>::setComponentReference
(
    const label patchi,
    const label facei,
    const direction,
    const scalar value
)
{
    if (psi_.needReference())
    {
        if (Pstream::master())
        {
            internalCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]];

            boundaryCoeffs_[patchi][facei] +=
                diag()[psi_.mesh().boundary()[patchi].faceCells()[facei]]
               *value;
        }
    }
}


template<>
Foam::lduSolverPerformance Foam::fvMatrix<Foam::scalar>::solve
(
    const dictionary& solverControls
)
{
    profilingTrigger profSolve("fvMatrix::solve_" + psi_.name());

    if (debug)
    {
        Info<< "fvMatrix<scalar>::solve(const dictionary& solverControls) : "
               "solving fvMatrix<scalar>"
            << endl;
    }

    // Complete matrix assembly.  HJ, 17/Apr/2012
    completeAssembly();

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    scalarField totalSource = source_;
    addBoundarySource(totalSource, false);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Cast into a non-const to solve.  HJ, 6/May/2016
    GeometricField<scalar, fvPatchField, volMesh>& psi =
        const_cast<GeometricField<scalar, fvPatchField, volMesh>&>(psi_);

    // Solver call
    lduSolverPerformance solverPerf = lduSolver::New
    (
        psi_.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces,
        solverControls
    )->solve(psi.internalField(), totalSource);

    solverPerf.print();

    // Diagonal has been restored, clear complete assembly flag?
    diag() = saveDiag;

    psi.correctBoundaryConditions();

    psi_.mesh().solutionDict().setSolverPerformance(psi_.name(), solverPerf);

    return solverPerf;
}


template<>
Foam::tmp<Foam::scalarField> Foam::fvMatrix<Foam::scalar>::residual() const
{
    // Complete matrix assembly.  HJ, 17/Apr/2012
    fvMatrix<scalar>& m = const_cast<fvMatrix<scalar>&>(*this);
    m.completeAssembly();

    // To avoid copying the diagonal, execute boundary diag multiplication
    // separately.  See source provided to lduMatrix::residual
    // HJ, 26/Jun/2015
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    // Bug fix: boundary source must be added on rhs
    // HJ, 26/Jun/2015
    scalarField totalSource = source_;
    addBoundarySource(totalSource, false);

   // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.internalField(),
            totalSource - boundaryDiag*psi_.internalField(),
            boundaryCoeffs_,
            interfaces,
            0
        )
    );

    return tres;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H() const
{
    tmp<volScalarField> tHphi
    (
        new volScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimVol,
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& Hphi = tHphi();

    Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().V();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


template<>
Foam::tmp<Foam::volScalarField> Foam::fvMatrix<Foam::scalar>::H1() const
{
    tmp<volScalarField> tH1
    (
        new volScalarField
        (
            IOobject
            (
                "H(1)",
                psi_.instance(),
                psi_.mesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/(dimVol*psi_.dimensions()),
            zeroGradientFvPatchScalarField::typeName
        )
    );
    volScalarField& H1_ = tH1();

    H1_.internalField() = lduMatrix::H1();
    //addBoundarySource(Hphi.internalField());

    H1_.internalField() /= psi_.mesh().V();
    H1_.correctBoundaryConditions();

    return tH1;
}


// ************************************************************************* //
