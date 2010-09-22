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
     Finite-Area scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "faScalarMatrix.H"
#include "zeroGradientFaPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Set reference level for a component of the solution
// on a given patch face
template<>
void faMatrix<scalar>::setComponentReference
(
    const label patchI,
    const label edgeI,
    const direction,
    const scalar value
)
{
    const unallocLabelList& faceLabels =
        psi_.mesh().boundary()[patchI].edgeFaces();

    internalCoeffs_[patchI][edgeI] +=
        diag()[faceLabels[edgeI]];

    boundaryCoeffs_[patchI][edgeI] = value;
}


template<>
autoPtr<faMatrix<scalar>::faSolver> faMatrix<scalar>::solver
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<scalar>::solver(const dictionary&) : "
               "solver for faMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    autoPtr<faMatrix<scalar>::faSolver> solverPtr
    (
        new faMatrix<scalar>::faSolver
        (
            *this,
            lduSolver::New
            (
                psi_.name(),
                *this,
                boundaryCoeffs_,
                internalCoeffs_,
                interfaces,
                solverControls
            )
        )
    );

    diag() = saveDiag;

    return solverPtr;
}


template<>
lduSolverPerformance faMatrix<scalar>::faSolver::solve
(
    const dictionary& solverControls
)
{
    scalarField saveDiag = faMat_.diag();
    faMat_.addBoundaryDiag(faMat_.diag(), 0);

    scalarField totalSource = faMat_.source();
    faMat_.addBoundarySource(totalSource, false);

    solver_->read(solverControls);
    lduSolverPerformance solverPerf =
        solver_->solve(faMat_.psi().internalField(), totalSource);

    solverPerf.print();

    faMat_.diag() = saveDiag;

    faMat_.psi().correctBoundaryConditions();

    return solverPerf;
}


template<>
lduSolverPerformance faMatrix<scalar>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "faMatrix<scalar>::solve(const dictionary&) : "
               "solving faMatrix<scalar>"
            << endl;
    }

    scalarField saveDiag = diag();
    addBoundaryDiag(diag(), 0);

    scalarField totalSource = source_;
    addBoundarySource(totalSource, 0);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Solver call
    lduSolverPerformance solverPerf = lduSolver::New
    (
        psi_.name(),
        *this,
        boundaryCoeffs_,
        internalCoeffs_,
        interfaces,
        solverControls
    )->solve(psi_.internalField(), totalSource);

    solverPerf.print();

    diag() = saveDiag;

    psi_.correctBoundaryConditions();

    return solverPerf;
}


// Return the matrix residual
template<>
tmp<scalarField> faMatrix<scalar>::residual() const
{
    scalarField boundaryDiag(psi_.size(), 0.0);
    addBoundaryDiag(boundaryDiag, 0);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.internalField(),
            source_ - boundaryDiag*psi_.internalField(),
            boundaryCoeffs_,
            interfaces,
            0
        )
    );

    addBoundarySource(tres());

    return tres;
}


// H operator
template<>
tmp<areaScalarField> faMatrix<scalar>::H() const
{
    tmp<areaScalarField> tHphi
    (
        new areaScalarField
        (
            IOobject
            (
                "H("+psi_.name()+')',
                psi_.instance(),
                psi_.db(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            psi_.mesh(),
            dimensions_/dimArea,
            zeroGradientFaPatchScalarField::typeName
        )
    );
    areaScalarField Hphi = tHphi();

    Hphi.internalField() = (lduMatrix::H(psi_.internalField()) + source_);
    addBoundarySource(Hphi.internalField());

    Hphi.internalField() /= psi_.mesh().S();
    Hphi.correctBoundaryConditions();

    return tHphi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
