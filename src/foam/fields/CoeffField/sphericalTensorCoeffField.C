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

#include "sphericalTensorCoeffField.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CoeffField<Foam::sphericalTensor>::CoeffField(const label size)
:
    DecoupledCoeffField<sphericalTensor>(size)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const CoeffField<sphericalTensor>& f
)
:
    DecoupledCoeffField<sphericalTensor>(f)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const DecoupledCoeffField<sphericalTensor>& f
)
:
    DecoupledCoeffField<sphericalTensor>(f)
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField
(
    const tmp<DecoupledCoeffField<sphericalTensor> >& tf
)
:
    DecoupledCoeffField<sphericalTensor>(tf())
{}


Foam::CoeffField<Foam::sphericalTensor>::CoeffField(Istream& is)
:
    DecoupledCoeffField<sphericalTensor>(is)
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor> >& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>::scalarTypeField& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor>::scalarTypeField>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const CoeffField<sphericalTensor>::linearTypeField& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


void Foam::CoeffField<Foam::sphericalTensor>::operator=
(
    const tmp<CoeffField<sphericalTensor>::linearTypeField>& f
)
{
    DecoupledCoeffField<sphericalTensor>::operator=(f);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const CoeffField<sphericalTensor>& f
)
{
    const DecoupledCoeffField<sphericalTensor>& df = f;
    return operator<<(os, df);
}


Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const tmp<CoeffField<sphericalTensor> >& tf
)
{
    const DecoupledCoeffField<sphericalTensor>& df = tf();
    os << df;
    tf.clear();
    return os;
}


/* * * * * * * * * * * * * * * * Global functions  * * * * * * * * * * * * * */

Foam::tmp<Foam::CoeffField<Foam::sphericalTensor> > Foam::inv
(
    const CoeffField<sphericalTensor>& f
)
{
    const DecoupledCoeffField<sphericalTensor>& df = f;

    return tmp<CoeffField<sphericalTensor> >
    (
        new CoeffField<sphericalTensor>(inv(df)())
    );
}


template<>
void Foam::multiply
(
    sphericalTensorField& f,
    const CoeffField<sphericalTensor>& f1,
    const sphericalTensor& f2
)
{
    const DecoupledCoeffField<sphericalTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    sphericalTensorField& f,
    const CoeffField<sphericalTensor>& f1,
    const sphericalTensorField& f2
)
{
    const DecoupledCoeffField<sphericalTensor>& df1 = f1;

    multiply(f, df1, f2);
}


template<>
void Foam::multiply
(
    sphericalTensorField& f,
    const sphericalTensorField& f1,
    const CoeffField<sphericalTensor>& f2
)
{
    const DecoupledCoeffField<sphericalTensor>& df2 = f2;

    multiply(f, f1, df2);
}


// ************************************************************************* //
