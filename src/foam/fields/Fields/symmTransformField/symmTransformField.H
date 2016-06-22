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

Typedef
    Foam::symmTransformField

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    symmTransformField.C

\*---------------------------------------------------------------------------*/

#ifndef symmTransformField_H
#define symmTransformField_H

#include "symmTransform.H"
#include "symmTensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void transform(Field<Type>&, const symmTensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const symmTensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const symmTensorField&, const tmp<Field<Type> >&);

template<class Type>
tmp<Field<Type> > transform(const tmp<symmTensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tmp<symmTensorField>&, const tmp<Field<Type> >&);


template<class Type>
void transform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const tmp<Field<Type> >&);


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const symmTensorField&
);

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<symmTensorField>&
);


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const symmTensorField&
);

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<symmTensorField>&
);


template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const symmTensorField&
);

template<>
tmp<Field<tensor> > transformFieldMask<tensor>
(
    const tmp<symmTensorField>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "symmTransformField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
