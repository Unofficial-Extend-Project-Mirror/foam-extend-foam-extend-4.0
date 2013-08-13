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
     Tet Finite Element scalar matrix member functions and operators

\*---------------------------------------------------------------------------*/

#include "tetFemScalarMatrix.H"
#include "tetPointFields.H"
#include "tetPolyPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<>
lduSolverPerformance tetFemMatrix<scalar>::solve
(
    const dictionary& solverControls
)
{
    if (debug)
    {
        Info<< "tetFemMatrix<scalar>::solve(const dictionary&) : "
            << "solving tetFemMatrix<scalar>"
            << endl;
    }

    // Add boundary source for gradient-type conditions
    addBoundarySourceDiag();

    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    // Set component boundary conditions
    scalarField sourceCpy = source_;

    // Store the boundary coefficients for insertion of boundary conditions
    setComponentBoundaryConditions(0, psi_, sourceCpy);

    // Add the coupling coefficients
    addCouplingCoeffs();
    addCouplingSource(sourceCpy);

    // prepare for coupled interface update
    FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());
    FieldField<Field, scalar> coupledIntCoeffs(psi_.boundaryField().size());

    forAll(psi_.boundaryField(), patchI)
    {
        const tetPolyPatchScalarField& ptf = psi_.boundaryField()[patchI];

        coupledBouCoeffs.set
        (
            patchI,
            new scalarField(ptf.cutBouCoeffs(*this))
        );

        coupledIntCoeffs.set
        (
            patchI,
            new scalarField(ptf.cutIntCoeffs(*this))
        );
    }

    eliminateCouplingCoeffs();

    lduSolverPerformance solverPerf
    (
        "tetFemMatrix<scalar>::solve",
        psi_.name()
    );

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    solverPerf = lduSolver::New
    (
        psi_.name(),
        *this,
        coupledBouCoeffs,
        coupledIntCoeffs,
        interfaces,
        solverControls
    )->solve(psi_.internalField(), sourceCpy);

    solverPerf.print();

    reconstructMatrix();

    if (debug)
    {
        Info<< "tetFemMatrix<scalar>::solve(const scalar, const scalar) : "
            << "correcting boundary conditions"
            << endl;
    }

    psi_.correctBoundaryConditions();

    return solverPerf;
}


// Return the matrix residual
template<>
tmp<scalarField> tetFemMatrix<scalar>::residual()
{
    // Store the boundary coefficients for insertion of boundary conditions
    storeBoundaryCoeffs();

    // Set component boundary conditions
    scalarField sourceCpy = source_;

    setComponentBoundaryConditions(0, psi_, sourceCpy);

    // Add the coupling coefficients
    addCouplingCoeffs();
    addCouplingSource(sourceCpy);

    // Prepare for coupled interface update
    FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());

    forAll(psi_.boundaryField(), patchI)
    {
        const tetPolyPatchScalarField& ptf = psi_.boundaryField()[patchI];
        coupledBouCoeffs.set
        (
            patchI,
            new scalarField(ptf.cutBouCoeffs(*this))
        );
    }

    eliminateCouplingCoeffs();

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    tmp<scalarField> tres
    (
        lduMatrix::residual
        (
            psi_.internalField(),
            sourceCpy,
            coupledBouCoeffs,
            interfaces,
            0
        )
    );

    reconstructMatrix();

    return tres;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
