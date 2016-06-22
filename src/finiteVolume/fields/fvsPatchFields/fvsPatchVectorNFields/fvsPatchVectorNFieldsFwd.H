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

#ifndef fvsPatchVectorNFieldsFwd_H
#define fvsPatchVectorNFieldsFwd_H

#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type> class fvsPatchField;

typedef fvsPatchField<vector2> fvsPatchVector2Field;
typedef fvsPatchField<vector3> fvsPatchVector3Field;
typedef fvsPatchField<vector4> fvsPatchVector4Field;
typedef fvsPatchField<vector6> fvsPatchVector6Field;
typedef fvsPatchField<vector8> fvsPatchVector8Field;

typedef fvsPatchField<tensor2> fvsPatchTensor2Field;
typedef fvsPatchField<tensor3> fvsPatchTensor3Field;
typedef fvsPatchField<tensor4> fvsPatchTensor4Field;
typedef fvsPatchField<tensor6> fvsPatchTensor6Field;
typedef fvsPatchField<tensor8> fvsPatchTensor8Field;

typedef fvsPatchField<diagTensor2> fvsPatchDiagTensor2Field;
typedef fvsPatchField<diagTensor3> fvsPatchDiagTensor3Field;
typedef fvsPatchField<diagTensor4> fvsPatchDiagTensor4Field;
typedef fvsPatchField<diagTensor6> fvsPatchDiagTensor6Field;
typedef fvsPatchField<diagTensor8> fvsPatchDiagTensor8Field;

typedef fvsPatchField<sphericalTensor2> fvsPatchSphericalTensor2Field;
typedef fvsPatchField<sphericalTensor3> fvsPatchSphericalTensor3Field;
typedef fvsPatchField<sphericalTensor4> fvsPatchSphericalTensor4Field;
typedef fvsPatchField<sphericalTensor6> fvsPatchSphericalTensor6Field;
typedef fvsPatchField<sphericalTensor8> fvsPatchSphericalTensor8Field;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
