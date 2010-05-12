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
    Multiply a given vector (second argument) by the matrix or its transpose
    and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::Amul
(
    scalarField& Ax,
    const scalarField& x,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    // Reset multiplication result to zero
    // HJ, 5/Nov/2007
    Ax = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        coupleBouCoeffs,
        interfaces,
        x,
        Ax,
        cmpt
    );

    // AmulCore must be additive to account for initialisation step
    // in ldu interfaces.  HJ, 6/Nov/2007
    AmulCore(Ax, x);

    // Update coupled interfaces
    updateMatrixInterfaces
    (
        coupleBouCoeffs,
        interfaces,
        x,
        Ax,
        cmpt
    );
}


void Foam::lduMatrix::AmulCore
(
    scalarField& Ax,
    const scalarField& x
) const
{
    scalar* __restrict__ AxPtr = Ax.begin();

    const scalar* const __restrict__ xPtr = x.begin();

    // Protection for multiplication of incomplete matrices
    // HJ, 19/Sep/2008
    if (hasDiag())
    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        register const label nCells = diag().size();
        for (register label cell=0; cell<nCells; cell++)
        {
            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&xPtr[cell+96],0,1);
            __builtin_prefetch (&diagPtr[cell+96],0,1);
            __builtin_prefetch (&AxPtr[cell+96],1,1);
            #endif

            // AmulCore must be additive to account for initialisation step
            // in ldu interfaces.  HJ, 6/Nov/2007
            AxPtr[cell] += diagPtr[cell]*xPtr[cell];
        }
    }

    // Protection for multiplication of incomplete matrices
    // HJ, 19/Sep/2008
    if (hasUpper() || hasLower())
    {
        const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* const __restrict__ upperPtr = upper().begin();
        const scalar* const __restrict__ lowerPtr = lower().begin();

        register const label nFaces = upper().size();
        #ifdef ICC_IA64_PREFETCH
        #pragma swp
        #endif
        for (register label face=0; face<nFaces; face++)
        {
            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&uPtr[face+32],0,0);
            __builtin_prefetch (&lPtr[face+32],0,0);
            __builtin_prefetch (&lowerPtr[face+32],0,1);
            __builtin_prefetch (&xPtr[lPtr[face+32]],0,1);
            __builtin_prefetch (&AxPtr[uPtr[face+32]],0,1);
            #endif

            AxPtr[uPtr[face]] += lowerPtr[face]*xPtr[lPtr[face]];

            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&upperPtr[face+32],0,1);
            __builtin_prefetch (&xPtr[uPtr[face+32]],0,1);
            __builtin_prefetch (&AxPtr[lPtr[face+32]],0,1);
            #endif

            AxPtr[lPtr[face]] += upperPtr[face]*xPtr[uPtr[face]];
        }
    }
}


void Foam::lduMatrix::Tmul
(
    scalarField& Tx,
    const scalarField& x,
    const FieldField<Field, scalar>& coupleIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    // Reset multiplication result to zero
    // HJ, 5/Nov/2007
    Tx = 0;

    // Initialise the update of coupled interfaces
    initMatrixInterfaces
    (
        coupleIntCoeffs,
        interfaces,
        x,
        Tx,
        cmpt
    );

    // TmulCore must be additive to account for initialisation step
    // in ldu interfaces.  HJ, 6/Nov/2007
    TmulCore(Tx, x);

    // Update coupled interfaces
    updateMatrixInterfaces
    (
        coupleIntCoeffs,
        interfaces,
        x,
        Tx,
        cmpt
    );
}


void Foam::lduMatrix::TmulCore
(
    scalarField& Tx,
    const scalarField& x
) const
{
    scalar* __restrict__ TxPtr = Tx.begin();

    const scalar* const __restrict__ xPtr = x.begin();

    // Protection for multiplication of incomplete matrices
    // HJ, 19/Sep/2008
    if (hasDiag())
    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        register const label nCells = diag().size();
        for (register label cell=0; cell<nCells; cell++)
        {
            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&xPtr[cell+96],0,1);
            __builtin_prefetch (&diagPtr[cell+96],0,1);
            __builtin_prefetch (&TxPtr[cell+96],1,1);
            #endif

            // TmulCore must be additive to account for initialisation step
            // in ldu interfaces.  HJ, 6/Nov/2007
            TxPtr[cell] += diagPtr[cell]*xPtr[cell];
        }
    }

    // Protection for multiplication of incomplete matrices
    // HJ, 19/Sep/2008
    if (hasUpper() || hasLower())
    {
        const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* const __restrict__ lowerPtr = lower().begin();
        const scalar* const __restrict__ upperPtr = upper().begin();

        register const label nFaces = upper().size();
        for (register label face=0; face<nFaces; face++)
        {
            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&uPtr[face+32],0,0);
            __builtin_prefetch (&lPtr[face+32],0,0);
            __builtin_prefetch (&upperPtr[face+32],0,1);
            __builtin_prefetch (&xPtr[lPtr[face+32]],0,1);
            __builtin_prefetch (&TxPtr[uPtr[face+32]],0,1);
            #endif

            TxPtr[uPtr[face]] += upperPtr[face]*xPtr[lPtr[face]];

            #ifdef ICC_IA64_PREFETCH
            __builtin_prefetch (&lowerPtr[face+32],0,1);
            __builtin_prefetch (&xPtr[uPtr[face+32]],0,1);
            __builtin_prefetch (&TxPtr[lPtr[face+32]],0,1);
            #endif

            TxPtr[lPtr[face]] += lowerPtr[face]*xPtr[uPtr[face]];
        }
    }
}


void Foam::lduMatrix::sumA
(
    scalarField& sumA,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces
) const
{
    scalar* __restrict__ sumAPtr = sumA.begin();

    const scalar* __restrict__ diagPtr = diag().begin();

    const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
    const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

    const scalar* __restrict__ lowerPtr = lower().begin();
    const scalar* __restrict__ upperPtr = upper().begin();

    register const label nCells = diag().size();
    register const label nFaces = upper().size();

    for (register label cell=0; cell<nCells; cell++)
    {
        #ifdef ICC_IA64_PREFETCH
        __builtin_prefetch (&diagPtr[cell+96],0,1);
        __builtin_prefetch (&sumAPtr[cell+96],1,1);
        #endif

        sumAPtr[cell] = diagPtr[cell];
    }

    #ifdef ICC_IA64_PREFETCH
    #pragma swp  
    #endif

    for (register label face=0; face<nFaces; face++)
    {
        #ifdef ICC_IA64_PREFETCH
        __builtin_prefetch (&uPtr[face+32],0,0);
        __builtin_prefetch (&lPtr[face+32],0,0);
        __builtin_prefetch (&lowerPtr[face+32],0,1);
        __builtin_prefetch (&sumAPtr[uPtr[face+32]],0,1);
        #endif

        sumAPtr[uPtr[face]] += lowerPtr[face];

        #ifdef ICC_IA64_PREFETCH
        __builtin_prefetch (&upperPtr[face+32],0,1);
        __builtin_prefetch (&sumAPtr[lPtr[face+32]],0,1);
        #endif

        sumAPtr[lPtr[face]] += upperPtr[face];
    }

    // Add the interface internal coefficients to diagonal
    // and the interface boundary coefficients to the sum-off-diagonal
    forAll(interfaces, patchI)
    {
        if (interfaces.set(patchI))
        {
            const unallocLabelList& pa = lduAddr().patchAddr(patchI);
            const scalarField& pCoeffs = coupleBouCoeffs[patchI];

            forAll(pa, face)
            {
                sumAPtr[pa[face]] -= pCoeffs[face];
            }
        }
    }
}


void Foam::lduMatrix::residual
(
    scalarField& rA,
    const scalarField& x,
    const scalarField& b,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    // Reset multiplication result to zero
    // HJ, 5/Nov/2007
    rA = 0;

    // Standard implementation
    Amul(rA, x, coupleBouCoeffs, interfaces, cmpt);

    const scalar* const __restrict__ bPtr = b.begin();
    scalar* __restrict__ rAPtr = rA.begin();

    forAll (rA, cell)
    {
        #ifdef ICC_IA64_PREFETCH
        __builtin_prefetch (&bPtr[cell+96],0,1);
        __builtin_prefetch (&rAPtr[cell+96],1,1);
        #endif

        rAPtr[cell] = bPtr[cell] - rAPtr[cell];
    }

    // Optimised access.  Obsolete, since AMul is fully optimised.
    // HJ, 5-6/Nov/2007
//     scalar* __restrict__ rAPtr = rA.begin();

//     const scalar* const __restrict__ xPtr = x.begin();
//     const scalar* const __restrict__ diagPtr = diag().begin();
//     const scalar* const __restrict__ bPtr = b.begin();

//     const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
//     const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

//     const scalar* const __restrict__ upperPtr = upper().begin();
//     const scalar* const __restrict__ lowerPtr = lower().begin();

//     // Parallel boundary initialisation.
//     // Note: there is a change of sign in the coupled
//     // interface update.  The reason for this is that the
//     // internal coefficients are all located at the l.h.s. of
//     // the matrix whereas the "implicit" coefficients on the
//     // coupled boundaries are all created as if the
//     // coefficient contribution is of a b-kind (i.e. they
//     // have a sign as if they are on the r.h.s. of the matrix.
//     // To compensate for this, it is necessary to turn the
//     // sign of the contribution.

//     FieldField<Field, scalar> mBouCoeffs(coupleBouCoeffs.size());

//     forAll(mBouCoeffs, patchi)
//     {
//         if (interfaces.set(patchi))
//         {
//             mBouCoeffs.set(patchi, -coupleBouCoeffs[patchi]);
//         }
//     }

//     // Initialise the update of coupled interfaces
//     initMatrixInterfaces
//     (
//         mBouCoeffs,
//         interfaces,
//         x,
//         rA, 
//         cmpt
//     );

//     register const label nCells = diag().size();
//     for (register label cell=0; cell<nCells; cell++)
//     {
//         #ifdef ICC_IA64_PREFETCH
//         __builtin_prefetch (&xPtr[cell+96],0,1);
//         __builtin_prefetch (&diagPtr[cell+96],0,1);
//         __builtin_prefetch (&bPtr[cell+96],0,1);
//         __builtin_prefetch (&rAPtr[cell+96],1,1);
//         #endif

//         rAPtr[cell] = bPtr[cell] - diagPtr[cell]*xPtr[cell];
//     }


//     register const label nFaces = upper().size();
//     #ifdef ICC_IA64_PREFETCH
//     #pragma swp  
//     #endif

//     for (register label face=0; face<nFaces; face++)
//     {
//         #ifdef ICC_IA64_PREFETCH
//         __builtin_prefetch (&uPtr[face+32],0,0);
//         __builtin_prefetch (&lPtr[face+32],0,0);
//         __builtin_prefetch (&lowerPtr[face+32],0,1);
//         __builtin_prefetch (&xPtr[lPtr[face+32]],0,1);
//         __builtin_prefetch (&rAPtr[uPtr[face+32]],0,1);
//         #endif

//         rAPtr[uPtr[face]] -= lowerPtr[face]*xPtr[lPtr[face]];

//         #ifdef ICC_IA64_PREFETCH
//         __builtin_prefetch (&upperPtr[face+32],0,1);
//         __builtin_prefetch (&xPtr[uPtr[face+32]],0,1);
//         __builtin_prefetch (&rAPtr[lPtr[face+32]],0,1);
//         #endif

//         rAPtr[lPtr[face]] -= upperPtr[face]*xPtr[uPtr[face]];
//     }

//     // Update coupled interfaces
//     updateMatrixInterfaces
//     (
//         mBouCoeffs,
//         interfaces,
//         x,
//         rA,
//         cmpt
//     );
//     Info << "res: " << gSumMag(rA) << endl;
}


Foam::tmp<Foam::scalarField> Foam::lduMatrix::residual
(
    const scalarField& x,
    const scalarField& b,
    const FieldField<Field, scalar>& coupleBouCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const direction cmpt
) const
{
    tmp<scalarField> trA(new scalarField(x.size()));
    residual(trA(), x, b, coupleBouCoeffs, interfaces, cmpt);
    return trA;
}


// ************************************************************************* //
