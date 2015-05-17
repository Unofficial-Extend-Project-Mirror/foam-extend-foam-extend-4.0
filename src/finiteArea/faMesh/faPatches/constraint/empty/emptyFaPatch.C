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

Description

\*---------------------------------------------------------------------------*/

#include "emptyFaPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Patch name
defineTypeNameAndDebug(emptyFaPatch, 0);

// Add the patch constructor functions to the hash tables
addToRunTimeSelectionTable(faPatch, emptyFaPatch, dictionary);

// Over-riding the face normals return from the underlying patch
// This is the only piece of info used out of the underlying primitivePatch
// I choose to store it there because it is used in primitive patch operations
// and it should not be duplicated as before.  However, to ensure everything
// in the empty patch is sized to zero, we shall here return a regerence to
// a zero-sized field (it does not matter what the field is
//
// const vectorField& emptyFaPatch::edgeNormals() const
// {
//     return faceAreas();
// }


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
