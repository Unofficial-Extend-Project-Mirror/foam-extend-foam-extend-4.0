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


// ************************************************************************* //
