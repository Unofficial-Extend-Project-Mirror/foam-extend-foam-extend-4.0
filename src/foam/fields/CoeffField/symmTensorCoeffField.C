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

#include "symmTensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::symmTensor>::CoeffField(const label size)
:
    DecoupledCoeffField<symmTensor>(size)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const CoeffField<symmTensor>& f
)
:
    DecoupledCoeffField<symmTensor>(f)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const DecoupledCoeffField<symmTensor>& f
)
:
    DecoupledCoeffField<symmTensor>(f)
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField
(
    const tmp<DecoupledCoeffField<symmTensor> >& tf
)
:
    DecoupledCoeffField<symmTensor>(tf())
{}


Foam::CoeffField<Foam::symmTensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<symmTensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor> >& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const CoeffField<symmTensor>::linearTypeField& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


void Foam::CoeffField<Foam::symmTensor>::operator=
(
    const tmp<CoeffField<symmTensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<symmTensor>::operator=(f);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CoeffField<symmTensor>& f)
{
    const DecoupledCoeffField<symmTensor>& df = f;
    return operator<<(os, df);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<symmTensor> >& tf
)
{
    const DecoupledCoeffField<symmTensor>& df = tf();
    os << df;
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::symmTensor> > Foam::inv
(
    const CoeffField<symmTensor>& f
)
{
    const DecoupledCoeffField<symmTensor>& df = f;

    return tmp<CoeffField<symmTensor> >
    (
        new CoeffField<symmTensor>(inv(df)())
    );
}


template<>
void Foam::multiply
(
    symmTensorField& f,
    const CoeffField<symmTensor>& f1,
    const symmTensor& f2
)
{
    const DecoupledCoeffField<symmTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    symmTensorField& f,
    const CoeffField<symmTensor>& f1,
    const symmTensorField& f2
)
{
    const DecoupledCoeffField<symmTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    symmTensorField& f,
    const symmTensorField& f1,
    const CoeffField<symmTensor>& f2
)
{
    const DecoupledCoeffField<symmTensor>& df2 = f2;

    multiply(f, f1, df2);
}


// ************************************************************************* //
