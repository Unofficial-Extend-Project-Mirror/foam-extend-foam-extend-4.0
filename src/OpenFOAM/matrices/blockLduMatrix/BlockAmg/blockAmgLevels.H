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

#include "BlockAmgLevel.H"
#include "fineBlockAmgLevel.H"
#include "coarseBlockAmgLevel.H"

// KRJ included to make the different types
#include "VectorTensorNFieldsFwd.H"
#include "ExpandTensorNField.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// KRJ: 2012-12-15: Is this enough to have the templating?

typedef BlockAmgLevel<scalar> blockAmgScalarLevel;
typedef BlockAmgLevel<vector> blockAmgVectorLevel;
typedef BlockAmgLevel<tensor> blockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef BlockAmgLevel<type > blockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

typedef coarseBlockAmgLevel<scalar> coarseblockAmgScalarLevel;
typedef coarseBlockAmgLevel<vector> coarseblockAmgVectorLevel;
typedef coarseBlockAmgLevel<tensor> coarseblockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef coarseBlockAmgLevel<type > coarseblockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

typedef fineBlockAmgLevel<scalar> fineblockAmgScalarLevel;
typedef fineBlockAmgLevel<vector> fineblockAmgVectorLevel;
typedef fineBlockAmgLevel<tensor> fineblockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef fineBlockAmgLevel<type > fineblockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
