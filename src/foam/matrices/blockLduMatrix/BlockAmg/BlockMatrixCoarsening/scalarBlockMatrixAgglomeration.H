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
    BlockMatrixAgglomeration

Description
    Specialisation of the BlockMatrixAgglomeration for scalars.

Author
    Klas Jareteg, 2013-01-31

\*---------------------------------------------------------------------------*/

#ifndef scalarBlockMatrixAgglomeration_H
#define scalarBlockMatrixAgglomeration_H

#include "blockMatrixCoarsenings.H"
#include "blockMatrixAgglomerations.H"
#include "BlockMatrixAgglomeration.H"
#include "BlockMatrixCoarsening.H"
#include "runTimeSelectionTables.H"
#include "scalarBlockLduMatrix.H"
#include "scalarBlockConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * Forward declaration of template friend fuctions * * * * * * * //


/*---------------------------------------------------------------------------*\
                        Class BlockMatrixAgglomeration Declaration
\*---------------------------------------------------------------------------*/

//- Restrict matrix

template<>
inline autoPtr<BlockLduMatrix<scalar> >
BlockMatrixAgglomeration<scalar>::restrictMatrix() const
{
    FatalErrorIn
    (
        "autoPtr<BlockLduMatrix<scalar> > "
        "BlockMatrixAgglomeration<Type>::restrictMatrix() const"
    )   << "Function not implemented for Type=scalar. " << endl
        << abort(FatalError);

    // Dummy return to keep compiler happy
    return autoPtr<BlockLduMatrix<scalar> >(NULL);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
