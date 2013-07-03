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

#include "blockAmgLevels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineNamedTemplateTypeNameAndDebug(blockAmgScalarLevel, 0);
defineNamedTemplateTypeNameAndDebug(blockAmgVectorLevel, 0);
defineNamedTemplateTypeNameAndDebug(blockAmgTensorLevel, 0);
#define makeNamedTemplateTypeNameAndDebug(type, Type, args...)      \
    defineNamedTemplateTypeNameAndDebug(blockAmg##Type##Level, 0);

    forAllVectorNTypes(makeNamedTemplateTypeNameAndDebug);

#undef makeNamedTemplateTypeNameAndDebug

defineNamedTemplateTypeNameAndDebug(coarseblockAmgScalarLevel, 0);
defineNamedTemplateTypeNameAndDebug(coarseblockAmgVectorLevel, 0);
defineNamedTemplateTypeNameAndDebug(coarseblockAmgTensorLevel, 0);
#define makeNamedTemplateTypeNameAndDebug(type, Type, args...)      \
    defineNamedTemplateTypeNameAndDebug(coarseblockAmg##Type##Level, 0);

    forAllVectorNTypes(makeNamedTemplateTypeNameAndDebug);

#undef makeNamedTemplateTypeNameAndDebug

defineNamedTemplateTypeNameAndDebug(fineblockAmgScalarLevel, 0);
defineNamedTemplateTypeNameAndDebug(fineblockAmgVectorLevel, 0);
defineNamedTemplateTypeNameAndDebug(fineblockAmgTensorLevel, 0);
#define makeNamedTemplateTypeNameAndDebug(type, Type, args...)      \
    defineNamedTemplateTypeNameAndDebug(fineblockAmg##Type##Level, 0);

    forAllVectorNTypes(makeNamedTemplateTypeNameAndDebug);

#undef makeNamedTemplateTypeNameAndDebug

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
