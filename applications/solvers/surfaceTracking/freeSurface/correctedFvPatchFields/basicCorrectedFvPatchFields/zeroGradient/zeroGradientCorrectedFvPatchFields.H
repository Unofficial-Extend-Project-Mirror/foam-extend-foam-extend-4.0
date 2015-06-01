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
    zeroGradientCorrectedFvPatchFields

Description

SourceFiles
    zeroGradientCorrectedFvPatchFields.C


\*---------------------------------------------------------------------------*/

#ifndef zeroGradientCorrectedFvPatchFields_H
#define zeroGradientCorrectedFvPatchFields_H

#include "zeroGradientCorrectedFvPatchField.H"
#include "fieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef zeroGradientCorrectedFvPatchField<scalar> zeroGradientCorrectedFvPatchScalarField;
typedef zeroGradientCorrectedFvPatchField<vector> zeroGradientCorrectedFvPatchVectorField;
typedef zeroGradientCorrectedFvPatchField<tensor> zeroGradientCorrectedFvPatchTensorField;
typedef zeroGradientCorrectedFvPatchField<sphericalTensor> zeroGradientCorrectedFvPatchSphericalTensorField;
typedef zeroGradientCorrectedFvPatchField<symmTensor> zeroGradientCorrectedFvPatchSymmTensorField;
typedef zeroGradientCorrectedFvPatchField<symmTensor4thOrder>
    zeroGradientCorrectedFvPatchSymmTensor4thOrderField;
typedef zeroGradientCorrectedFvPatchField<diagTensor>
    zeroGradientCorrectedFvPatchDiagTensorField;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
