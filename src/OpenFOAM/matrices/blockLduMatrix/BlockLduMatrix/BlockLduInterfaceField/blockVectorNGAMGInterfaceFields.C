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
