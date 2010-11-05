/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
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
typedef fvsPatchField<vector4> fvsPatchVector4Field;
typedef fvsPatchField<vector6> fvsPatchVector6Field;
typedef fvsPatchField<vector8> fvsPatchVector8Field;

typedef fvsPatchField<tensor2> fvsPatchTensor2Field;
typedef fvsPatchField<tensor4> fvsPatchTensor4Field;
typedef fvsPatchField<tensor6> fvsPatchTensor6Field;
typedef fvsPatchField<tensor8> fvsPatchTensor8Field;

typedef fvsPatchField<diagTensor2> fvsPatchDiagTensor2Field;
typedef fvsPatchField<diagTensor4> fvsPatchDiagTensor4Field;
typedef fvsPatchField<diagTensor6> fvsPatchDiagTensor6Field;
typedef fvsPatchField<diagTensor8> fvsPatchDiagTensor8Field;

typedef fvsPatchField<sphericalTensor2> fvsPatchSphericalTensor2Field;
typedef fvsPatchField<sphericalTensor4> fvsPatchSphericalTensor4Field;
typedef fvsPatchField<sphericalTensor6> fvsPatchSphericalTensor6Field;
typedef fvsPatchField<sphericalTensor8> fvsPatchSphericalTensor8Field;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
