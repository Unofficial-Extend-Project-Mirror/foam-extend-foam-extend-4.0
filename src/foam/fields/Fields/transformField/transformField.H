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

InNamespace
    Foam

Description
    Spatial transformation functions for primitive fields.

SourceFiles
    transformField.C

\*---------------------------------------------------------------------------*/

#ifndef transformField_H
#define transformField_H

#include "transform.H"
#include "quaternion.H"
#include "septernion.H"
#include "vectorField.H"
#include "tensorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void transform(Field<Type>&, const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensorField&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensorField&, const tmp<Field<Type> >&);

template<class Type>
tmp<Field<Type> > transform(const tmp<tensorField>&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tmp<tensorField>&, const tmp<Field<Type> >&);


template<class Type>
void transform(Field<Type>&, const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const Field<Type>&);

template<class Type>
tmp<Field<Type> > transform(const tensor&, const tmp<Field<Type> >&);


template<class Type1, class Type2>
tmp<Field<Type1> > transformFieldMask(const Field<Type2>&);

template<class Type1, class Type2>
tmp<Field<Type1> > transformFieldMask(const tmp<Field<Type2> >&);


template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tensorField&
);

template<>
tmp<Field<symmTensor> > transformFieldMask<symmTensor>
(
    const tmp<tensorField>&
);


template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tensorField&
);

template<>
tmp<Field<sphericalTensor> > transformFieldMask<sphericalTensor>
(
    const tmp<tensorField>&
);

template<>
tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
(
    const tensorField&
 );

template<>
tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
(
    const tmp<tensorField>&
 );

template<>
tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
(
     const symmTensorField&
 );

 template<>
 tmp<Field<symmTensor4thOrder> > transformFieldMask<symmTensor4thOrder>
 (
    const tmp<symmTensorField>&
  );

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
     const tensorField&
 );

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
    const tmp<tensorField>&
 );

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
     const symmTensorField&
 );

template<>
tmp<Field<diagTensor> > transformFieldMask<diagTensor>
(
     const tmp<symmTensorField>&
 );


//- Rotate given vectorField with the given quaternion
void transform(vectorField&, const quaternion&, const vectorField&);

//- Rotate given vectorField with the given quaternion
tmp<vectorField> transform(const quaternion&, const vectorField&);

//- Rotate given tmp<vectorField> with the given quaternion
tmp<vectorField> transform(const quaternion&, const tmp<vectorField>&);


//- Transform given vectorField with the given septernion
void transform(vectorField&, const septernion&, const vectorField&);

//- Transform given vectorField with the given septernion
tmp<vectorField> transform(const septernion&, const vectorField&);

//- Transform given tmp<vectorField> with the given septernion
tmp<vectorField> transform(const septernion&, const tmp<vectorField>&);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "transformFieldTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
