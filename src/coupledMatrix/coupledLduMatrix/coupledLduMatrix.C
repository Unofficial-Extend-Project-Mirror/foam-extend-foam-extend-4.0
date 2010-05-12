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

Class
    coupledLduMatrix

Description
    Collection of lduMatrices solved together as a block system

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "coupledLduMatrix.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(coupledLduMatrix, 1);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given size
Foam::coupledLduMatrix::coupledLduMatrix(const label size)
:
    PtrList<lduMatrix>(size)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::coupledLduMatrix::~coupledLduMatrix()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::coupledLduMatrix::diagonal() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool diag = true;

    forAll (matrices, matrixI)
    {
        diag = diag && matrices[matrixI].diagonal();
    }

    return diag;
}


bool Foam::coupledLduMatrix::symmetric() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool sym = true;

    forAll (matrices, matrixI)
    {
        sym =
            (sym && matrices[matrixI].diagonal())
         || (sym && matrices[matrixI].symmetric());
    }

    return sym;
}


bool Foam::coupledLduMatrix::asymmetric() const
{
    const PtrList<lduMatrix>& matrices = *this;

    bool asym = false;

    forAll (matrices, matrixI)
    {
        asym = (asym || matrices[matrixI].asymmetric());
    }

    return asym;
}


void Foam::coupledLduMatrix::Amul
(
    FieldField<Field, scalar>& result,
    const FieldField<Field, scalar>& x,
    const PtrList<FieldField<Field, scalar> >& bouCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    // Reset product to zero
    result = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        bouCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );

    forAll (matrices, rowI)
    {
        matrices[rowI].AmulCore(result[rowI], x[rowI]);
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        bouCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );
}


void Foam::coupledLduMatrix::Tmul
(
    FieldField<Field, scalar>& result,
    const FieldField<Field, scalar>& x,
    const PtrList<FieldField<Field, scalar> >& intCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    // Reset product to zero
    result = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        intCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );

    forAll (matrices, rowI)
    {
        matrices[rowI].TmulCore(result[rowI], x[rowI]);
    }

    // Update couple interfaces
    updateMatrixInterfaces
    (
        intCoeffs,
        interfaces,
        x,
        result,
        cmpt
    );
}


void Foam::coupledLduMatrix::initMatrixInterfaces
(
    const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const FieldField<Field, scalar>& x,
    FieldField<Field, scalar>& result,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    forAll (matrices, rowI)
    {
        matrices[rowI].initMatrixInterfaces
        (
            coupleCoeffs[rowI],
            interfaces[rowI],
            x[rowI],
            result[rowI],
            cmpt
        );
    }
}


void Foam::coupledLduMatrix::updateMatrixInterfaces
(
    const PtrList<FieldField<Field, scalar> >& coupleCoeffs,
    const lduInterfaceFieldPtrsListList& interfaces,
    const FieldField<Field, scalar>& x,
    FieldField<Field, scalar>& result,
    const direction cmpt
) const
{
    const PtrList<lduMatrix>& matrices = *this;

    forAll (matrices, rowI)
    {
        matrices[rowI].updateMatrixInterfaces
        (
            coupleCoeffs[rowI],
            interfaces[rowI],
            x[rowI],
            result[rowI],
            cmpt
        );
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
