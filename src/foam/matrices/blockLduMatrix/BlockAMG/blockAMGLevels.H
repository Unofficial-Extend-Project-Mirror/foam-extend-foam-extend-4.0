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

#include "BlockAMGLevel.H"
#include "fineBlockAMGLevel.H"
#include "coarseBlockAMGLevel.H"

// KRJ included to make the different types
#include "VectorTensorNFieldsFwd.H"
#include "ExpandTensorNField.H"
#include "VectorNFieldTypes.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

// KRJ: 2012-12-15: Is this enough to have the templating?

typedef BlockAMGLevel<scalar> blockAmgScalarLevel;
typedef BlockAMGLevel<vector> blockAmgVectorLevel;
typedef BlockAMGLevel<tensor> blockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef BlockAMGLevel<type > blockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

typedef coarseBlockAMGLevel<scalar> coarseblockAmgScalarLevel;
typedef coarseBlockAMGLevel<vector> coarseblockAmgVectorLevel;
typedef coarseBlockAMGLevel<tensor> coarseblockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef coarseBlockAMGLevel<type > coarseblockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

typedef fineBlockAMGLevel<scalar> fineblockAmgScalarLevel;
typedef fineBlockAMGLevel<vector> fineblockAmgVectorLevel;
typedef fineBlockAMGLevel<tensor> fineblockAmgTensorLevel;

#define makeTypedef(type, Type, args...)                \
    typedef fineBlockAMGLevel<type > fineblockAmg##Type##Level;
forAllVectorNTypes(makeTypedef)
#undef makeTypedef

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
