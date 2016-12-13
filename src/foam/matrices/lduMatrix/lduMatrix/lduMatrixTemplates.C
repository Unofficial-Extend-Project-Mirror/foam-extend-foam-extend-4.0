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
    lduMatrix member H operations.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::lduMatrix::H(const Field<Type>& psi) const
{
    tmp<Field<Type> > tHpsi
    (
        new Field<Type>(lduAddr().size(), pTraits<Type>::zero)
    );

    if (lowerPtr_ || upperPtr_)
    {
        Field<Type> & Hpsi = tHpsi();

        Type* __restrict__ HpsiPtr = Hpsi.begin();

        const Type* __restrict__ psiPtr = psi.begin();

        const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
        const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

        const scalar* __restrict__ lowerPtr = lower().begin();
        const scalar* __restrict__ upperPtr = upper().begin();

        register const label nFaces = upper().size();

        for (register label face=0; face<nFaces; face++)
        {
            HpsiPtr[uPtr[face]] -= lowerPtr[face]*psiPtr[lPtr[face]];
            HpsiPtr[lPtr[face]] -= upperPtr[face]*psiPtr[uPtr[face]];
        }
    }

    return tHpsi;
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::lduMatrix::H(const tmp<Field<Type> >& tpsi) const
{
    tmp<Field<Type> > tHpsi(H(tpsi()));
    tpsi.clear();
    return tHpsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::lduMatrix::faceH(const Field<Type>& psi) const
{
    tmp<Field<Type> > tfaceHpsi
    (
        new Field<Type> (lduAddr().lowerAddr().size())
    );
    Field<Type>& faceHpsi = tfaceHpsi();

    if (lowerPtr_ || upperPtr_)
    {
        const scalarField& Lower = const_cast<const lduMatrix&>(*this).lower();
        const scalarField& Upper = const_cast<const lduMatrix&>(*this).upper();

        const unallocLabelList& l = lduAddr().lowerAddr();
        const unallocLabelList& u = lduAddr().upperAddr();

        for (register label face=0; face<l.size(); face++)
        {
            faceHpsi[face] = Upper[face]*psi[u[face]]
                - Lower[face]*psi[l[face]];
        }
    }
    else
    {
        // No off-diagonal.  Bug fix for conjugate matrices.  HJ, 27/Nov/2008
        faceHpsi = pTraits<Type>::zero;
    }

    return tfaceHpsi;
}


template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::lduMatrix::faceH(const tmp<Field<Type> >& tpsi) const
{
    tmp<Field<Type> > tfaceHpsi(faceH(tpsi()));
    tpsi.clear();
    return tfaceHpsi;
}


// ************************************************************************* //
