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

Class
    starMeshReader

Description
    Read pro-STAR vrt/cel/bnd files.

 * Starting with Pro-STAR version 4, the files have become easier to read.
 * Vertices are space-delimited.
 * The cell format is logical.
 * Trimmed and degenerate cells are saved as polyhedral.
 * The boundaries corresponds to cells and their faces.

 The protected data in meshReader are filled.

SourceFiles
    starMeshReader.C

\*---------------------------------------------------------------------------*/

#ifndef starMeshReader_H
#define starMeshReader_H

#include "meshReader.H"
#include "cellShape.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class starMeshReader Declaration
\*---------------------------------------------------------------------------*/

class starMeshReader : public meshReader
{
    // Private data

    //- Points supporting the mesh
    pointField points_;

    //- Cell shapes
    cellShapeList cellShapes_;

    //- Point labels (imported Point numbering not necessarily contiguous)
    labelList mapToFoamPointId_;

    //- Cell labels (imported Cell numbering not necessarily contiguous)
    labelList mapToFoamCellId_;

    // Private Member Functions

    //- Disallow default bitwise copy construct
    starMeshReader(const starMeshReader&);

    //- Disallow default bitwise assignment
    void operator=(const starMeshReader&);

public:

    // Constructors

    //- Construct from case name
    starMeshReader
    (
        const fileName& prefix,
        const Time& runtime,
        const scalar scaleFactor = 1.0
    );

    // Destructor
    virtual ~starMeshReader();

    // Member Functions
    //
    //- Return mesh points
    virtual const pointField & points() const
    {
        return points_;
    }

private:

    //- Read the mesh from the file(s)
    virtual void readMesh(const fileName &prefix, const scalar scaleFactor);
    void readPoints(const fileName &, const scalar scaleFactor);
    void readCells(const fileName &);
    void readBoundary(const fileName &);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
