/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    BlockILUC0Precon

Description
    Typedefs for ILUCp preconditioning

Author
    Vuko Vukcevic, FMENA Zagreb. All rights reserved.

SourceFiles
    blockILUC0Precons.C

\*---------------------------------------------------------------------------*/

#ifndef blockILUC0Precons_H
#define blockILUC0Precons_H

// Specialisations not implemented. VV, 3/Jul/2015.
#include "scalarBlockILUC0Precon.H"
#include "tensorBlockILUC0Precon.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef BlockILUC0Precon<scalar> blockILUC0PreconScalar;
typedef BlockILUC0Precon<vector> blockILUC0PreconVector;
typedef BlockILUC0Precon<tensor> blockILUC0PreconTensor;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
