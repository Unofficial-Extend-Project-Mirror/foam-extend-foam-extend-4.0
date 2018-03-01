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
            const overlapFringe& overlapFringeAlgorithm,
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
