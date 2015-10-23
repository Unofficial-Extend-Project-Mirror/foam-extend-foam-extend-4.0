/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Tetrahedral Finite Element matrix check for diagonal dominance.
    Especially useful for debugging complex interfaces

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class Type>
void tetFemMatrix<Type>::check()
{
    if (debug)
    {
        Info<< "tetFemMatrix<Type>::check : checking tetFemMatrix<Type>"
            << endl;
    }

    Pout << "First check diagonal dominance" << endl;

    // Get constant matrix access
    const tetFemMatrix<Type>& constMatrix = *this;

    const scalarField& lower = constMatrix.lower();
    const scalarField& upper = constMatrix.upper();
    const scalarField& dd = constMatrix.diag();

    // Calculate local matrix off-diag sum
    scalarField dummySource(lduAddr().size(), 0.0);

    const scalarField oldLower = lower;
    const scalarField oldUpper = upper;
    const scalarField oldDiag = dd;

    // Get matrix addressing
    const unallocLabelList& L = lduAddr().lowerAddr();
    const unallocLabelList& U = lduAddr().upperAddr();

    {
        scalarField matrixSumOffDiag(lduAddr().size(), 0.0);
        scalarField dummyInternal(lduAddr().size(), 1.0);

        forAll (L, face)
        {
            matrixSumOffDiag[L[face]] += oldLower[face];
            matrixSumOffDiag[U[face]] += oldUpper[face];
        }

        Pout<< "void tetFemMatrix<Type>::check() : "
            << "Raw matrix difference: "
            << sum(mag(matrixSumOffDiag + oldDiag)) << endl;
    }

    // Add the coupling coefficients
    addCouplingCoeffs();

    addCouplingSource(dummySource);

    // Make a copy of interfaces: no longer a reference
    // HJ, 20/Nov/2007
    lduInterfaceFieldPtrsList interfaces = psi_.boundaryField().interfaces();

    // Prepare for coupled interface update
    FieldField<Field, scalar> coupledBouCoeffs(psi_.boundaryField().size());
    FieldField<Field, scalar> coupledIntCoeffs(psi_.boundaryField().size());

    forAll(psi_.boundaryField(), patchI)
    {
        const tetPolyPatchField<Type>& ptf = psi_.boundaryField()[patchI];

        coupledBouCoeffs.set(patchI, ptf.cutBouCoeffs(*this));

        coupledIntCoeffs.set(patchI, ptf.cutIntCoeffs(*this));
    }

    eliminateCouplingCoeffs();

    Pout << "Second check diagonal dominance" << endl;
    // Calculate local matrix off-diag sum
    {
        scalarField matrixSumOffDiag(lduAddr().size(), 0.0);
        scalarField dummyInternal(lduAddr().size(), 1.0);

        forAll (L, face)
        {
            matrixSumOffDiag[L[face]] += lower[face];
            matrixSumOffDiag[U[face]] += upper[face];
        }

        // Initialise coupling
        forAll (interfaces, interfaceI)
        {
            // Check if the pointer is set
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].initInterfaceMatrixUpdate
                (
                    dummyInternal,
                    matrixSumOffDiag,
                    *this,
                    coupledBouCoeffs[interfaceI],
                    0,
                    static_cast<Pstream::commsTypes>
                    (
                        Pstream::defaultCommsType()
                    ),
                    false                       // Do not switch to lhs
                );
            }
        }

        // Update coupled interface
        forAll (interfaces, interfaceI)
        {
            if (interfaces.set(interfaceI))
            {
                interfaces[interfaceI].updateInterfaceMatrix
                (
                    dummyInternal,
                    matrixSumOffDiag,
                    *this,
                    coupledBouCoeffs[interfaceI],
                    0,
                    static_cast<Pstream::commsTypes>(Pstream::defaultCommsType()),
                    false                       // Do not switch to lhs
                );
            }
        }

        Pout<< "void tetFemMatrix<Type>::check() : "
            << "parallel matrix difference: "
            << sum(mag(matrixSumOffDiag + dd)) << endl;
    }

    // restore the matrix in the original state
    if (hasUpper())
    {
        this->upper() = oldUpper;
    }

    if (hasLower())
    {
        this->lower() = oldLower;
    }

    if (hasDiag())
    {
        this->diag() = oldDiag;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
