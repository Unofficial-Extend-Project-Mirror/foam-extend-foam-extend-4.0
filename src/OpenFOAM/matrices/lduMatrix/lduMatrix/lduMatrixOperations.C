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
    lduMatrix member operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::sumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    for (register label odcI = 0; odcI < l.size(); odcI++)
    {
        Diag[l[odcI]] += Lower[odcI];
        Diag[u[odcI]] += Upper[odcI];
    }
}


void Foam::lduMatrix::negSumDiag()
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();
    scalarField& Diag = diag();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    for (register label odcI = 0; odcI < l.size(); odcI++)
    {
        Diag[l[odcI]] -= Lower[odcI];
        Diag[u[odcI]] -= Upper[odcI];
    }
}


void Foam::lduMatrix::sumMagOffDiag
(
    scalarField& sumOff
) const
{
    const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
    const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

    const unallocLabelList& l = lduAddr().lowerAddr();
    const unallocLabelList& u = lduAddr().upperAddr();

    for (register label odcI = 0; odcI < l.size(); odcI++)
    {
        sumOff[u[odcI]] += mag(Lower[odcI]);
        sumOff[l[odcI]] += mag(Upper[odcI]);
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::operator=(const lduMatrix& A)
{
    if (this == &A)
    {
        FatalError
            << "lduMatrix::operator=(const lduMatrix&) : "
            << "attempted assignment to self"
            << abort(FatalError);
    }

    if (A.lowerPtr_)
    {
        lower() = A.lower();
    }
    else if (lowerPtr_)
    {
        delete lowerPtr_;
        lowerPtr_ = NULL;
    }

    if (A.upperPtr_)
    {
        upper() = A.upper();
    }
    else if (upperPtr_)
    {
        delete upperPtr_;
        upperPtr_ = NULL;
    }

    if (A.diagPtr_)
    {
        diag() = A.diag();
    }
}


void Foam::lduMatrix::negate()
{
    if (lowerPtr_)
    {
        lowerPtr_->negate();
    }

    if (upperPtr_)
    {
        upperPtr_->negate();
    }

    if (diagPtr_)
    {
        diagPtr_->negate();
    }
}


void Foam::lduMatrix::operator+=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() += A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() += A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() += A.upper();
        lower() += A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() += A.upper();
            upper() += A.upper();
        }
        else
        {
            lower() += A.lower();
            upper() += A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() += A.lower();
        upper() += A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorIn("lduMatrix::operator+=(const lduMatrix& A)")
            << "Unknown matrix type combination"
            << abort(FatalError);
    }
}


void Foam::lduMatrix::operator-=(const lduMatrix& A)
{
    if (A.diagPtr_)
    {
        diag() -= A.diag();
    }

    if (symmetric() && A.symmetric())
    {
        upper() -= A.upper();
    }
    else if (symmetric() && A.asymmetric())
    {
        if (upperPtr_)
        {
            lower();
        }
        else
        {
            upper();
        }

        upper() -= A.upper();
        lower() -= A.lower();
    }
    else if (asymmetric() && A.symmetric())
    {
        if (A.upperPtr_)
        {
            lower() -= A.upper();
            upper() -= A.upper();
        }
        else
        {
            lower() -= A.lower();
            upper() -= A.lower();
        }

    }
    else if (asymmetric() && A.asymmetric())
    {
        lower() -= A.lower();
        upper() -= A.upper();
    }
    else if (diagonal())
    {
        if (A.upperPtr_)
        {
            upper() = -A.upper();
        }

        if (A.lowerPtr_)
        {
            lower() = -A.lower();
        }
    }
    else if (A.diagonal())
    {
    }
    else
    {
        FatalErrorIn("lduMatrix::operator-=(const lduMatrix& A)")
            << "Unknown matrix type combination"
            << abort(FatalError);
    }
}


void Foam::lduMatrix::operator*=(const scalarField& sf)
{
    if (diagPtr_)
    {
        *diagPtr_ *= sf;
    }

    if (upperPtr_)
    {
        scalarField& upper = *upperPtr_;

        const unallocLabelList& l = lduAddr().lowerAddr();

        for (register label odcI = 0; odcI < upper.size(); odcI++)
        {
            upper[odcI] *= sf[l[odcI]];
        }
    }

    if (lowerPtr_)
    {
        scalarField& lower = *lowerPtr_;

        const unallocLabelList& u = lduAddr().upperAddr();

        for (register label odcI = 0; odcI < lower.size(); odcI++)
        {
            lower[odcI] *= sf[u[odcI]];
        }
    }
}


void Foam::lduMatrix::operator*=(scalar s)
{
    if (diagPtr_)
    {
        *diagPtr_ *= s;
    }

    if (upperPtr_)
    {
        *upperPtr_ *= s;
    }

    if (lowerPtr_)
    {
        *lowerPtr_ *= s;
    }
}


Foam::tmp<Foam::scalarField > Foam::lduMatrix::H1() const
{
    tmp<scalarField > tH1
    (
        new scalarField(lduAddr().size(), 0.0)
    );

    if (lowerPtr_ || upperPtr_)
    {
        scalarField& H1_ = tH1();

        scalar* __restrict__ H1Ptr = H1_.begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        register const label nOdcIs = upper().size();

        for (register label odcI = 0; odcI < nOdcIs; odcI++)
        {
            H1Ptr[uPtr[odcI]] -= lowerPtr[odcI];
            H1Ptr[lPtr[odcI]] -= upperPtr[odcI];
        }
    }

    return tH1;
}


// ************************************************************************* //
