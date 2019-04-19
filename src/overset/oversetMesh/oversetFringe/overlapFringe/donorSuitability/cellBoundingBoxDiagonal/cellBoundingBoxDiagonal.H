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
    Foam::donorSuitability::cellBoundingBoxDiagonal

Description
    Cell bounding box diagonal suitability. Donor suitability function
    corresponds to the diagonal of a bounding box encompassing a cell.

    Minimises the overlap by shifting it towards cells with similar bounding box
    diagonal (example: a tet cell with vertices on a hex cell faces (opposite
    faces in a hex, have the same bounding box). This might be suitable for
    unstructured tet meshes.

SourceFiles
    cellBoundingBoxDiagonal.C

Author
    Vuko Vukcevic, FMENA Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef cellBoundingBoxDiagonal_H
#define cellBoundingBoxDiagonal_H

#include "donorSuitability.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace donorSuitability
{

/*---------------------------------------------------------------------------*\
                    Class cellBoundingBoxDiagonal Declaration
\*---------------------------------------------------------------------------*/

class cellBoundingBoxDiagonal
:
    public donorSuitability
{

public:

    //- Runtime type information
    TypeName("cellBoundingBoxDiagonal");


    // Constructors

        //- Construct from components
        cellBoundingBoxDiagonal
        (
            const oversetFringe& oversetFringeAlgorithm,
            const dictionary& dict
        );


    //- Destructor
    virtual ~cellBoundingBoxDiagonal()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace cellBoundingBoxDiagonal
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
