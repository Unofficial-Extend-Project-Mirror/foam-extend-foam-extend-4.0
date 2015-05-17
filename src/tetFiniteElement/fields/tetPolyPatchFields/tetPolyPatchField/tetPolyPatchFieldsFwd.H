/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

Class
    tetPolyPatchField

Description

\*---------------------------------------------------------------------------*/

#ifndef patchFieldsFwd_H
#define patchFieldsFwd_H

#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type> class tetPolyPatchField;

typedef tetPolyPatchField<scalar> tetPolyPatchScalarField;
typedef tetPolyPatchField<vector> tetPolyPatchVectorField;
typedef tetPolyPatchField<sphericalTensor> tetPolyPatchSphericalTensorField;
typedef tetPolyPatchField<symmTensor> tetPolyPatchSymmTensorField;
typedef tetPolyPatchField<symmTensor4thOrder> tetPolyPatchSymmTensor4thOrderField;
typedef tetPolyPatchField<diagTensor> tetPolyPatchDiagTensorField;
typedef tetPolyPatchField<tensor> tetPolyPatchTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
