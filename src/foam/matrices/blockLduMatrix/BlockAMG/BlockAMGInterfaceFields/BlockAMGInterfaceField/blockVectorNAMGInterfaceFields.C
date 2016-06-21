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

Class
    blockVectorNAMGInterfaceFields

Description
    Macros for VectorN types for AMG interface fields with block coeffs

Author
    Klas Jareteg, 2013-02-08

\*---------------------------------------------------------------------------*/

#include "BlockAMGInterfaceField.H"
#include "ProcessorBlockAMGInterfaceField.H"
#include "GGIBlockAMGInterfaceField.H"
#include "VectorNFieldTypes.H"
#include "ExpandTensorNField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
typedef BlockAMGInterfaceField<type > block##Type##AMGInterfaceField;         \
defineNamedTemplateTypeNameAndDebug(block##Type##AMGInterfaceField, 0);       \
defineTemplateRunTimeSelectionTable(block##Type##AMGInterfaceField, lduInterface); \
                                                                              \
typedef ProcessorBlockAMGInterfaceField<type > block##Type##ProcessorAMGInterfaceField;  \
makeBlockAMGInterfaceField(block##Type##AMGInterfaceField, block##Type##ProcessorAMGInterfaceField); \
typedef GGIBlockAMGInterfaceField<type > block##Type##GGIAMGInterfaceField;  \
makeBlockAMGInterfaceField(block##Type##AMGInterfaceField, block##Type##GGIAMGInterfaceField); \

forAllVectorNTypes(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
