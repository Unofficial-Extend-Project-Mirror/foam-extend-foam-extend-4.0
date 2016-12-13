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

Description

\*---------------------------------------------------------------------------*/

#include "AverageIOField.H"
#include "fieldTypes.H"

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

typedef AverageIOField<scalar> scalarAverageIOField;
typedef AverageIOField<vector> vectorAverageIOField;
typedef AverageIOField<sphericalTensor> sphericalTensorAverageIOField;
typedef AverageIOField<symmTensor> symmTensorAverageIOField;
typedef AverageIOField<symmTensor4thOrder> symmTensor4thOrderAverageIOField;
typedef AverageIOField<diagTensor> diagTensorAverageIOField;
typedef AverageIOField<tensor> tensorAverageIOField;

defineTemplateTypeNameAndDebugWithName
(
    scalarAverageIOField,
    "scalarAverageField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    vectorAverageIOField,
    "vectorAverageField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    sphericalTensorAverageIOField,
    "sphericalTensorAverageField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
    symmTensorAverageIOField,
    "symmTensorAverageField",
    0
);
defineTemplateTypeNameAndDebugWithName
(
 symmTensor4thOrderAverageIOField,
 "symmTensor4thOrderAverageField",
    0
 );
defineTemplateTypeNameAndDebugWithName
(
 diagTensorAverageIOField,
 "diagTensorAverageField",
    0
 );
defineTemplateTypeNameAndDebugWithName
(
    tensorAverageIOField,
    "tensorAverageField",
    0
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
