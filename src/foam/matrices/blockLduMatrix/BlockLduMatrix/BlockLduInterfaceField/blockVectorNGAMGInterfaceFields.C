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
    blockVectorNGAMGInterfaceFields

Description
    Macros for VectorN types for GAMG interface fields with block coeffs

Author
    Klas Jareteg, 2013-02-08

\*----------------------------------------------------------------------------*/

#include "BlockGAMGInterfaceField.H"
#include "processorBlockGAMGInterfaceField.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                                               \
                                                                                                        \
typedef BlockGAMGInterfaceField<type > block##Type##GAMGInterfaceField;                                 \
defineNamedTemplateTypeNameAndDebug(block##Type##GAMGInterfaceField, 0);                                \
defineTemplateRunTimeSelectionTable(block##Type##GAMGInterfaceField, lduInterface);                     \
                                                                                                        \
typedef processorBlockGAMGInterfaceField<type > block##Type##ProcessorGAMGInterfaceField;               \
makeBlockGAMGInterfaceField(block##Type##GAMGInterfaceField, block##Type##ProcessorGAMGInterfaceField); \

forAllVectorNTypes(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
