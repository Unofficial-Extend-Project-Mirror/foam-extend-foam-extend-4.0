/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    scalarCoeffField

Description

\*---------------------------------------------------------------------------*/

#include "scalarCoeffField.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

const char* const Foam::CoeffField<Foam::scalar>::typeName("CoeffField");


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::blockCoeffBase::activeLevel
Foam::CoeffField<Foam::scalar>::activeType() const
{
    return blockCoeffBase::SCALAR;
}


Foam::tmp<Foam::CoeffField<Foam::scalar> >
Foam::CoeffField<Foam::scalar>::transpose() const
{
    return tmp<CoeffField<scalar> >(new CoeffField<scalar>(*this));
}


const Foam::scalarField&
Foam::CoeffField<Foam::scalar>::asScalar() const
{
    return *this;
}


Foam::scalarField&
Foam::CoeffField<Foam::scalar>::asScalar()
{
    return *this;
}


const Foam::scalarField&
Foam::CoeffField<Foam::scalar>::asLinear() const
{
    return *this;
}


Foam::scalarField&
Foam::CoeffField<Foam::scalar>::asLinear()
{
    return *this;
}


Foam::BlockCoeff<Foam::scalar>
Foam::CoeffField<Foam::scalar>::getCoeff(const label index) const
{
    BlockCoeff<scalar> result;

    result.asScalar() = (operator[](index));

    return result;
}


void Foam::CoeffField<Foam::scalar>::setCoeff
(
    const label index,
    const BlockCoeff<scalar>& coeff
)
{
    operator[](index) = coeff.asScalar();
}


void Foam::CoeffField<Foam::scalar>::getSubset
(
    CoeffField<scalar>& f,
    const label start,
    const label size
) const
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorIn
        (
            "void Foam::CoeffField<Foam::scalar>::getSubset\n"
            "(\n"
            "    CoeffField<scalar>& f,\n"
            "    const label start,\n"
            "    const label size\n"
            ") const"
        )   << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    const scalarField& localF = *this;

    forAll (f, fI)
    {
        f[fI] = localF[start + fI];
    }
}


void Foam::CoeffField<Foam::scalar>::getSubset
(
    CoeffField<scalar>& f,
    const labelList& addr
) const
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "void Foam::CoeffField<Foam::scalar>::getSubset\n"
            "(\n"
            "    CoeffField<scalar>& f,\n"
            "    const labelList addr\n"
            ") const"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    const scalarField& localF = *this;

    forAll (f, fI)
    {
        f[fI] = localF[addr[fI]];
    }
}


void Foam::CoeffField<Foam::scalar>::setSubset
(
    const CoeffField<scalar>& f,
    const label start,
    const label size
)
{
    // Check sizes
    if (f.size() != size)
    {
        FatalErrorIn
        (
            "void Foam::CoeffField<Foam::scalar>::setSubset\n"
            "(\n"
            "     const CoeffField<scalar>& f,\n"
            "    const label start,\n"
            "    const label size\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << size
            << abort(FatalError);
    }

    scalarField& localF = *this;

    forAll (f, fI)
    {
        localF[start + fI] = f[fI];
    }
}


void Foam::CoeffField<Foam::scalar>::setSubset
(
    const CoeffField<scalar>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "void Foam::CoeffField<Foam::scalar>::setSubset\n"
            "(\n"
            "    const CoeffField<scalar>& f,\n"
            "    const labelList addr\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    scalarField& localF = this->asScalar();

    forAll (f, fI)
    {
        localF[addr[fI]] = f[fI];
    }
}


void Foam::CoeffField<Foam::scalar>::zeroOutSubset
(
    const label start,
    const label size
)
{
    scalarField& localF = *this;

    for (label ffI = 0; ffI < size; ffI++)
    {
        localF[start + ffI] = pTraits<scalar>::zero;
    }
}


void Foam::CoeffField<Foam::scalar>::zeroOutSubset
(
    const labelList& addr
)
{
    scalarField& localF = *this;

    forAll (addr, ffI)
    {
        localF[addr[ffI]] = pTraits<scalar>::zero;
    }
}


void Foam::CoeffField<Foam::scalar>::addSubset
(
    const CoeffField<scalar>& f,
    const labelList& addr
)
{
    // Check sizes
    if (f.size() != addr.size())
    {
        FatalErrorIn
        (
            "void Foam::CoeffField<Foam::scalar>::addSubset\n"
            "(\n"
            "    const CoeffField<scalar>& f,\n"
            "    const labelList addr\n"
            ")"
        )   << "Incompatible sizes: " << f.size() << " and " << addr.size()
            << abort(FatalError);
    }

    scalarField& localF = this->asScalar();

    forAll (f, fI)
    {
        localF[addr[fI]] += f[fI];
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::scalar>::operator=(const scalarField& f)
{
    scalarField::operator=(f);
}


void Foam::CoeffField<Foam::scalar>::operator=(const tmp<scalarField>& tf)
{
    scalarField::operator=(tf);
}


// * * * * * * * * * * * * * * * Ostream Operator  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CoeffField<scalar>& f)
{
    os << static_cast<const scalarField&>(f);

    return os;
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<scalar> >& tf
)
{
    os << tf();
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

template<>
Foam::tmp<Foam::CoeffField<Foam::scalar> >
Foam::inv(const CoeffField<scalar>& f)
{
    tmp<CoeffField<scalar> > tresult(new CoeffField<scalar>(f.size()));
    scalarField& result = tresult();
    result = 1.0/f;

    return tresult;
}


template<>
void Foam::negate
(
    CoeffField<scalar>& f,
    const CoeffField<scalar>& f1
)
{
    f = f1;
    f.negate();
}


template<>
void Foam::multiply
(
    scalarField& f,
    const CoeffField<scalar>& f1,
    const scalar& f2
)
{
    const scalarField& sf = f1;
    f = sf*f2;
}


template<>
void Foam::multiply
(
    scalarField& f,
    const CoeffField<scalar>& f1,
    const scalarField& f2
)
{
    const scalarField& sf = f1;
    f = sf*f2;
}


template<>
void Foam::multiply
(
    scalarField& f,
    const scalarField& f1,
    const CoeffField<scalar>& f2
)
{
    const scalarField& sf = f2;
    f = f1*sf;
}


// ************************************************************************* //
