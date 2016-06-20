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

#include "symmTensor4thOrderCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::symmTensor4thOrder>::CoeffField(const label size)
:
    DecoupledCoeffField<symmTensor4thOrder>(size)
{}


Foam::CoeffField<Foam::symmTensor4thOrder>::CoeffField
(
    const CoeffField<symmTensor4thOrder>& f
)
:
    DecoupledCoeffField<symmTensor4thOrder>(f)
{}


Foam::CoeffField<Foam::symmTensor4thOrder>::CoeffField
(
    const DecoupledCoeffField<symmTensor4thOrder>& f
)
:
    DecoupledCoeffField<symmTensor4thOrder>(f)
{}


Foam::CoeffField<Foam::symmTensor4thOrder>::CoeffField
(
    const tmp<DecoupledCoeffField<symmTensor4thOrder> >& tf
)
:
    DecoupledCoeffField<symmTensor4thOrder>(tf())
{}


Foam::CoeffField<Foam::symmTensor4thOrder>::CoeffField(Istream& is)
:
    DecoupledCoeffField<symmTensor4thOrder>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const CoeffField<symmTensor4thOrder>& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const tmp<CoeffField<symmTensor4thOrder> >& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const CoeffField<symmTensor4thOrder>::scalarTypeField& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const tmp<CoeffField<symmTensor4thOrder>::scalarTypeField>& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const CoeffField<symmTensor4thOrder>::linearTypeField& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor4thOrder>::operator=
(
    const tmp<CoeffField<symmTensor4thOrder>::linearTypeField>& f
)
{
    DecoupledCoeffField<symmTensor4thOrder>::operator=(f);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CoeffField<symmTensor4thOrder>& f
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df = f;
    return operator<<(os, df);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<symmTensor4thOrder> >& tf
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df = tf();
    os << df;
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::symmTensor4thOrder> > Foam::inv
(
    const CoeffField<symmTensor4thOrder>& f
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df = f;

    return tmp<CoeffField<symmTensor4thOrder> >
    (
        new CoeffField<symmTensor4thOrder>(inv(df)())
    );
}


template<>
void Foam::multiply
(
    symmTensor4thOrderField& f,
    const CoeffField<symmTensor4thOrder>& f1,
    const symmTensor4thOrder& f2
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    symmTensor4thOrderField& f,
    const CoeffField<symmTensor4thOrder>& f1,
    const symmTensor4thOrderField& f2
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    symmTensor4thOrderField& f,
    const symmTensor4thOrderField& f1,
    const CoeffField<symmTensor4thOrder>& f2
)
{
    const DecoupledCoeffField<symmTensor4thOrder>& df2 = f2;

    multiply(f, f1, df2);
}


// ************************************************************************* //
