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

#include "diagTensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::diagTensor>::CoeffField(const label size)
:
    DecoupledCoeffField<diagTensor>(size)
{}


Foam::CoeffField<Foam::diagTensor>::CoeffField
(
    const CoeffField<diagTensor>& f
)
:
    DecoupledCoeffField<diagTensor>(f)
{}


Foam::CoeffField<Foam::diagTensor>::CoeffField
(
    const DecoupledCoeffField<diagTensor>& f
)
:
    DecoupledCoeffField<diagTensor>(f)
{}


Foam::CoeffField<Foam::diagTensor>::CoeffField
(
    const tmp<DecoupledCoeffField<diagTensor> >& tf
)
:
    DecoupledCoeffField<diagTensor>(tf())
{}


Foam::CoeffField<Foam::diagTensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<diagTensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const CoeffField<diagTensor>& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const tmp<CoeffField<diagTensor> >& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const CoeffField<diagTensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const tmp<CoeffField<diagTensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const CoeffField<diagTensor>::linearTypeField& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


void Foam::CoeffField<Foam::diagTensor>::operator=
(
    const tmp<CoeffField<diagTensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<diagTensor>::operator=(f);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const CoeffField<diagTensor>& f)
{
    const DecoupledCoeffField<diagTensor>& df = f;
    return operator<<(os, df);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<diagTensor> >& tf
)
{
    const DecoupledCoeffField<diagTensor>& df = tf();
    os << df;
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::diagTensor> > Foam::inv
(
    const CoeffField<diagTensor>& f
)
{
    const DecoupledCoeffField<diagTensor>& df = f;

    return tmp<CoeffField<diagTensor> >
    (
        new CoeffField<diagTensor>(inv(df)())
    );
}


template<>
void Foam::multiply
(
    diagTensorField& f,
    const CoeffField<diagTensor>& f1,
    const diagTensor& f2
)
{
    const DecoupledCoeffField<diagTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    diagTensorField& f,
    const CoeffField<diagTensor>& f1,
    const diagTensorField& f2
)
{
    const DecoupledCoeffField<diagTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    diagTensorField& f,
    const diagTensorField& f1,
    const CoeffField<diagTensor>& f2
)
{
    const DecoupledCoeffField<diagTensor>& df2 = f2;

    multiply(f, f1, df2);
}


// ************************************************************************* //
