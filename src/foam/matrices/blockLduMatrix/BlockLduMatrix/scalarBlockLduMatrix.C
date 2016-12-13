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

#include "coeffFields.H"
#include "BlockLduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<>
void BlockLduMatrix<scalar>::sumDiag()
{
    scalarField& activeDiag = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] += activeUpper[coeffI];
            activeDiag[u[coeffI]] += activeUpper[coeffI];
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = lower();
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] += activeLower[coeffI];
            activeDiag[u[coeffI]] += activeUpper[coeffI];
        }
    }
    else
    {
        FatalErrorIn("void BlockLduMatrix<scalar>::sumDiag()")
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<>
void BlockLduMatrix<scalar>::negSumDiag()
{
    scalarField& activeDiag = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] -= activeUpper[coeffI];
            activeDiag[u[coeffI]] -= activeUpper[coeffI];
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = lower();
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiag[l[coeffI]] -= activeLower[coeffI];
            activeDiag[u[coeffI]] -= activeUpper[coeffI];
        }
    }
    else
    {
        FatalErrorIn("void BlockLduMatrix<scalar>::negSumDiag()")
            << "No off-diagonal available"
            << abort(FatalError);
    }
}


template<>
void BlockLduMatrix<scalar>::check() const
{
    // Copy the diagonal
    scalarField activeDiagCopy = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    if (symmetric())
    {
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiagCopy[l[coeffI]] -= activeUpper[coeffI];
            activeDiagCopy[u[coeffI]] -= activeUpper[coeffI];
        }

        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Symmetric matrix: raw matrix difference: "
            << sum(mag(activeDiagCopy))
            << " scaled: "
            << sum(mag(activeDiagCopy))/sum(mag(diag()))
            << endl;
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = lower();
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeDiagCopy[l[coeffI]] -= activeLower[coeffI];
            activeDiagCopy[u[coeffI]] -= activeUpper[coeffI];
        }

        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Asymmetric matrix: raw matrix difference: "
            << sum(mag(activeDiagCopy))
            << " scaled: "
            << sum(mag(activeDiagCopy))/sum(mag(diag()))
            << endl;
    }
    else
    {
        Info<< "void BlockLduMatrix<scalar>::check() const : "
            << "Diagonal matrix" << endl;
    }
}


template<>
void BlockLduMatrix<scalar>::relax
(
    const scalarField& x,
    scalarField& b,
    const scalar alpha
)
{
    scalarField& activeDiag = diag();

    scalarField activeDiagOld = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    scalarField sumOff(activeDiag.size(), 0.0);

    if (symmetric())
    {
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumOff[u[coeffI]] += mag(activeUpper[coeffI]);
            sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
        }
    }
    else if (asymmetric())
    {
        const scalarField& activeLower = lower();
        const scalarField& activeUpper = upper();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            sumOff[u[coeffI]] += mag(activeLower[coeffI]);
            sumOff[l[coeffI]] += mag(activeUpper[coeffI]);
        }
    }

    activeDiag = max(activeDiag, sumOff);
    activeDiag *= 1.0/alpha;

    // Add the relaxation contribution to b
    b += (activeDiag - activeDiagOld)*x;
}


template<>
void BlockLduMatrix<scalar>::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        scalarField& activeUpper = *upperPtr_;

        const unallocLabelList& l = lduAddr().lowerAddr();

        for (register label coeffI = 0; coeffI < l.size(); coeffI++)
        {
            activeUpper[coeffI] *= sf[l[coeffI]];
        }
    }

    if (lowerPtr_)
    {
        scalarField& activeLower = *lowerPtr_;

        const unallocLabelList& u = lduAddr().upperAddr();

        for (register label coeffI = 0; coeffI < u.size(); coeffI++)
        {
            activeLower[coeffI] *= sf[u[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::AmulCore
(
    scalarField& Ax,
    const scalarField& x
) const
{
    // Note: pointer looping
    const label* const __restrict__ U = lduAddr().upperAddr().begin();
    const label* const __restrict__ L = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ X = x.begin();
    scalar* __restrict__ AX = Ax.begin();

    if (thereIsDiag())
    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        register const label nCells = diag().size();
        for (register label cell=0; cell<nCells; cell++)
        {
            // AmulCore must be additive to account for initialisation step
            // in ldu interfaces.  HJ, 6/Nov/2007
            AX[cell] += diagPtr[cell]*X[cell];
        }
    }

    if (symmetric())
    {
        if (thereIsUpper())
        {
            const scalar* const __restrict__ Upper = upper().begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                AX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
                AX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
            }
        }
        else
        {
            const scalar* const __restrict__ Lower = lower().begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                AX[U[coeffI]] += Lower[coeffI]*X[L[coeffI]];
                AX[L[coeffI]] += Lower[coeffI]*X[U[coeffI]];
            }
        }
    }
    else if (asymmetric())
    {
        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__  Upper = upper().begin();

        for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            AX[U[coeffI]] += Lower[coeffI]*X[L[coeffI]];
            AX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::TmulCore
(
    scalarField& Tx,
    const scalarField& x
) const
{
    // Note: pointer looping
    const label* const __restrict__ U = lduAddr().upperAddr().begin();
    const label* const __restrict__ L = lduAddr().lowerAddr().begin();

    const scalar* const __restrict__ X = x.begin();
    scalar* __restrict__ TX = Tx.begin();

    if (thereIsDiag())
    {
        const scalar* const __restrict__ diagPtr = diag().begin();

        register const label nCells = diag().size();
        for (register label cell=0; cell<nCells; cell++)
        {
            // AmulCore must be additive to account for initialisation step
            // in ldu interfaces.  HJ, 6/Nov/2007
            TX[cell] += diagPtr[cell]*X[cell];
        }
    }

    if (symmetric())
    {
        if (thereIsUpper())
        {
            const scalar* const __restrict__ Upper = upper().begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                TX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
                TX[L[coeffI]] += Upper[coeffI]*X[U[coeffI]];
            }
        }
        else
        {
            const scalar* const __restrict__ Lower = lower().begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                TX[U[coeffI]] += Lower[coeffI]*X[L[coeffI]];
                TX[L[coeffI]] += Lower[coeffI]*X[U[coeffI]];
            }
        }
    }
    else if (asymmetric())
    {
        const scalar* const __restrict__ Lower = lower().begin();
        const scalar* const __restrict__ Upper = upper().begin();

        for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
        {
            TX[U[coeffI]] += Upper[coeffI]*X[L[coeffI]];
            TX[L[coeffI]] += Lower[coeffI]*X[U[coeffI]];
        }
    }
}


template<>
void BlockLduMatrix<scalar>::segregateB
(
    scalarField&,
    const scalarField&
) const
{
    FatalErrorIn
    (
        "void BlockLduMatrix<scalar>::segregateB\n"
        "(\n"
        "    scalarField&,\n"
        "    const scalarField&\n"
        ") const"
    )   << "Requested decoupling of scalar matrix - never coupled"
        << abort(FatalError);
}


template<>
tmp<scalarField> BlockLduMatrix<scalar>::H(const scalarField& x) const
{
    tmp<scalarField> tresult(new scalarField(lduAddr().size(), 0));
    scalarField& result = tresult();

    if (thereIsUpper() || thereIsLower())
    {
        // Note: pointer looping
        const label* const __restrict__ U = lduAddr().upperAddr().begin();
        const label* const __restrict__ L = lduAddr().lowerAddr().begin();

        const scalar* const __restrict__ X = x.begin();

        if (symmetric())
        {
            if (thereIsUpper())
            {
                const scalar* const __restrict__ Upper = upper().begin();

                scalar* __restrict__ R = result.begin();

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[U[coeffI]] -= Upper[coeffI]*X[U[coeffI]];
                }

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[L[coeffI]] -= Upper[coeffI]*X[L[coeffI]];
                }
            }
            else if (thereIsLower())
            {
                const scalar* const __restrict__ Lower = lower().begin();

                scalar* __restrict__ R = result.begin();

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[U[coeffI]] -= Lower[coeffI]*X[U[coeffI]];
                }

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[L[coeffI]] -= Lower[coeffI]*X[L[coeffI]];
                }
            }
        }
        else
        {
            const scalar* const __restrict__ Lower = lower().begin();
            const scalar* const __restrict__ Upper = upper().begin();

            scalar* __restrict__ R = result.begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                R[U[coeffI]] -= Upper[coeffI]*X[U[coeffI]];
            }

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                R[L[coeffI]] -= Lower[coeffI]*X[L[coeffI]];
            }
        }
    }

    return tresult;
}


template<>
tmp<scalarField> BlockLduMatrix<scalar>::faceH(const scalarField& x) const
{
    tmp<scalarField> tresult(new scalarField(lduAddr().upperAddr().size(), 0));
    scalarField& result = tresult();

    if (thereIsUpper() || thereIsLower())
    {
        if (symmetric())
        {
            if (thereIsUpper())
            {
                // Note: pointer looping
                const label* const __restrict__ U = lduAddr().upperAddr().begin();
                const label* const __restrict__ L = lduAddr().lowerAddr().begin();

                const scalar* const __restrict__ Upper = upper().begin();
                const scalar* const __restrict__ X = x.begin();

                scalar* __restrict__ R = result.begin();

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[coeffI] = Upper[coeffI]*(X[U[coeffI]] - X[L[coeffI]]);
                }
            }
            else if (thereIsLower())
            {
                // Note: pointer looping
                const label* const __restrict__ U = lduAddr().upperAddr().begin();
                const label* const __restrict__ L = lduAddr().lowerAddr().begin();

                const scalar* const __restrict__ Lower = lower().begin();
                const scalar* const __restrict__ X = x.begin();

                scalar* __restrict__ R = result.begin();

                for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
                {
                    R[coeffI] = Lower[coeffI]*(X[U[coeffI]] - X[L[coeffI]]);
                }
            }
        }
        else
        {
            // Note: pointer looping
            const label* const __restrict__ U = lduAddr().upperAddr().begin();
            const label* const __restrict__ L = lduAddr().lowerAddr().begin();

            const scalar* const __restrict__ Lower = lower().begin();
            const scalar* const __restrict__ Upper = upper().begin();
            const scalar* const __restrict__ X = x.begin();

            scalar* __restrict__ R = result.begin();

            for (register label coeffI = 0; coeffI < upper().size(); coeffI++)
            {
                R[coeffI] = Upper[coeffI]*X[U[coeffI]] - Lower[coeffI]*X[L[coeffI]];
            }
        }
    }

    return tresult;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
