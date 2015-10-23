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
    Foam::thresholdCellFaces

Description
    Selects the mesh cell faces specified by a threshold value.
    Non-triangulated by default.

SourceFiles
    thresholdCellFaces.C

\*---------------------------------------------------------------------------*/

#ifndef thresholdCellFaces_H
#define thresholdCellFaces_H

#include "MeshedSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class thresholdCellFaces Declaration
\*---------------------------------------------------------------------------*/

class thresholdCellFaces
:
    public MeshedSurface<face>
{
    //- Private typedefs for convenience
        typedef MeshedSurface<face> MeshStorage;

        //- Reference to mesh
        const polyMesh& mesh_;

        //- For every face the original cell in mesh
        labelList meshCells_;

    // Private Member Functions

        void calculate
        (
            const scalarField&,
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate
        );

public:

    //- Runtime type information
    TypeName("thresholdCellFaces");


    // Constructors

        //- Construct from mesh, field and threshold value
        thresholdCellFaces
        (
            const polyMesh&,
            const scalarField&,
            const scalar lowerThreshold,
            const scalar upperThreshold,
            const bool triangulate = false
        );


    // Member Functions

        //- For every face original cell in mesh
        labelList& meshCells()
        {
            return meshCells_;
        }

        //- For every face original cell in mesh
        const labelList& meshCells() const
        {
            return meshCells_;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
