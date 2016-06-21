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
    BlockMatrixAgglomeration

Description
    Specialisation of the BlockMatrixAgglomeration for tensors.

Author
    Klas Jareteg, 2013-01-31

\*---------------------------------------------------------------------------*/

#ifndef tensorBlockMatrixAgglomeration_H
#define tensorBlockMatrixAgglomeration_H

#include "blockMatrixCoarsenings.H"
#include "blockMatrixAgglomerations.H"
#include "BlockMatrixAgglomeration.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class BlockMatrixAgglomeration Declaration
\*---------------------------------------------------------------------------*/

//- Restrict matrix

template<>
inline autoPtr<BlockAMGLevel<tensor> >
BlockMatrixAgglomeration<tensor>::restrictMatrix
() const
{
    FatalErrorIn
    (
        "autoPtr<BlockAMGLevel<tensor> > "
        "BlockMatrixAgglomeration<Type>::"
        "restrictMatrix() const"
    )   << "Function not implemented for Type=tensor. " << endl
        << abort(FatalError);

    // Dummy return to keep compiler happy
    return autoPtr<BlockAMGLevel<tensor> >(NULL);
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
