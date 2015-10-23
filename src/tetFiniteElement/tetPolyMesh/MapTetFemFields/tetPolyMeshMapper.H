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
    tetPolyMeshMapper

Description
    Class holds all the necessary information for mapping fields associated
    with tetPolyMesh.
SourceFiles
    tetPolyMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyMeshMapper_H
#define tetPolyMeshMapper_H

#include "tetPointMapper.H"
#include "pointMapper.H"
#include "faceMapper.H"
#include "cellMapper.H"
#include "tetPolyBoundaryMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMesh;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                 Class tetPolyMeshMapper Declaration
\*---------------------------------------------------------------------------*/

class tetPolyMeshMapper
{
    // Private data

        //- Reference to mesh
        const tetPolyMesh& mesh_;

        //- Point mapper
        pointMapper pointMap_;

        //- Face mapper
        faceMapper faceMap_;

        //- Cell mapper
        cellMapper cellMap_;

        //- Point mapper (needs to be re-done to insert faces and cells)
        tetPointMapper tetPointMap_;

        //- Boundary mapper
        tetPolyBoundaryMapper boundaryMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyMeshMapper(const tetPolyMeshMapper&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyMeshMapper&);


public:

    // Constructors

        //- Construct from components
        tetPolyMeshMapper
        (
            const tetPolyMesh& mesh,
            const mapPolyMesh& meshMap
        )
        :
            mesh_(mesh),
            pointMap_(meshMap),
            faceMap_(meshMap),
            cellMap_(meshMap),
            tetPointMap_(mesh, meshMap, pointMap_, faceMap_, cellMap_),
            boundaryMap_(mesh, meshMap, pointMap_, faceMap_)
        {}


    // Member Functions

        //- Return mesh
        const tetPolyMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to objectRegistry storing fields
        //  Needed in mappers.  HJ, 16/Oct/2008
        const objectRegistry& thisDb() const
        {
            return mesh_.thisDb();
        }

        //- Return point mapper
        const morphFieldMapper& pointMap() const
        {
            return tetPointMap_;
        }

        //- Return element mapper
        const morphFieldMapper& elementMap() const
        {
            return cellMap_;
        }

        //- Return boundary mapper
        const tetPolyBoundaryMapper& boundaryMap() const
        {
            return boundaryMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
