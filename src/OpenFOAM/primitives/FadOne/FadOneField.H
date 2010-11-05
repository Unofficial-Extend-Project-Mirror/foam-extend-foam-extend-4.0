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
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    FadOneField

Description
    Specialisation of Field<T> for FadOne

SourceFiles
    FadOneFields.C

\*---------------------------------------------------------------------------*/

#ifndef FadOneField_H
#define FadOneField_H

#include "FadOne.H"
#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<int nVars>
tmp<Field<FadOne<nVars> > > ValueFadOneField(const UList<scalar>&);


template<int nVars>
tmp<scalarField> FadOneValue(const Field<FadOne<nVars> >&);

template<int nVars>
void FadOneSetValue
(
    Field<FadOne<nVars> >&,
    const scalarField&
);

template<int nVars>
tmp<scalarField> FadOneDeriv
(
    const Field<FadOne<nVars> >&,
    const direction
);

template<int nVars>
void FadOneSetDeriv
(
    Field<FadOne<nVars> >&,
    const direction,
    const scalarField&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#ifdef NoRepository
#   include "FadOneField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
