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
    blockVectorNLduInterfaceFields

Description
    Instantiation of block ldu interfaces for VectorN types

\*---------------------------------------------------------------------------*/

#include "BlockLduInterfaceField.H"
#include "GGIBlockLduInterfaceField.H"
#include "MixingPlaneBlockLduInterfaceField.H"
#include "OverlapGGIBlockLduInterfaceField.H"
#include "ProcessorBlockLduInterfaceField.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

#define makeTemplateTypeNameAndDebug(type, Type, args...)                     \
                                                                              \
defineTemplateTypeNameAndDebug(BlockLduInterfaceField<type>, 0);              \
                                                                              \
defineTemplateTypeNameAndDebug(GGIBlockLduInterfaceField<type>, 0);           \
defineTemplateTypeNameAndDebug(MixingPlaneBlockLduInterfaceField<type>, 0);   \
defineTemplateTypeNameAndDebug(OverlapGGIBlockLduInterfaceField<type>, 0);    \
defineTemplateTypeNameAndDebug(ProcessorBlockLduInterfaceField<type>, 0);

forAllVectorNTypes(makeTemplateTypeNameAndDebug);

#undef makeTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
