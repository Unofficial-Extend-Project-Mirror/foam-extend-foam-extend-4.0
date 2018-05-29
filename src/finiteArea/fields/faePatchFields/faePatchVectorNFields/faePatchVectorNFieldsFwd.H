/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#ifndef faePatchVectorNFieldsFwd_H
#define faePatchVectorNFieldsFwd_H

#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type> class faePatchField;

typedef faePatchField<vector2> faePatchVector2Field;
typedef faePatchField<vector3> faePatchVector3Field;
typedef faePatchField<vector4> faePatchVector4Field;
typedef faePatchField<vector6> faePatchVector6Field;
typedef faePatchField<vector8> faePatchVector8Field;

typedef faePatchField<tensor2> faePatchTensor2Field;
typedef faePatchField<tensor3> faePatchTensor3Field;
typedef faePatchField<tensor4> faePatchTensor4Field;
typedef faePatchField<tensor6> faePatchTensor6Field;
typedef faePatchField<tensor8> faePatchTensor8Field;

typedef faePatchField<diagTensor2> faePatchDiagTensor2Field;
typedef faePatchField<diagTensor3> faePatchDiagTensor3Field;
typedef faePatchField<diagTensor4> faePatchDiagTensor4Field;
typedef faePatchField<diagTensor6> faePatchDiagTensor6Field;
typedef faePatchField<diagTensor8> faePatchDiagTensor8Field;

typedef faePatchField<sphericalTensor2> faePatchSphericalTensor2Field;
typedef faePatchField<sphericalTensor3> faePatchSphericalTensor3Field;
typedef faePatchField<sphericalTensor4> faePatchSphericalTensor4Field;
typedef faePatchField<sphericalTensor6> faePatchSphericalTensor6Field;
typedef faePatchField<sphericalTensor8> faePatchSphericalTensor8Field;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
