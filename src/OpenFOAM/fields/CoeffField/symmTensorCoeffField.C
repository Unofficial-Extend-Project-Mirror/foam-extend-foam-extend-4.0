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


// ************************************************************************* //
