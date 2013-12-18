/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    BlockAamgPolicy

Description
    Typedefs for block AMG policies.

Author
    Klas Jareteg, 2012-12-15

SourceFiles
    blockAamgPolicies.C

\*---------------------------------------------------------------------------*/

#ifndef blockAamgPolicies_H
#define blockAamgPolicies_H

// Disabled instantiation for primitive types by force. HJ, 15/Apr/2013
#include "scalarBlockAamgPolicy.H"
#include "vectorBlockAamgPolicy.H"
#include "tensorBlockAamgPolicy.H"
#include "BlockAamgPolicy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockAamgPolicy<scalar> blockAamgPolicyScalar;
typedef BlockAamgPolicy<vector> blockAamgPolicyVector;
typedef BlockAamgPolicy<tensor> blockAamgPolicyTensor;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
