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
    tetPolyMeshMapperFaceDecomp

Description
    Class holds all the necessary information for mapping fields associated
    with tetPolyMeshFaceDecomp.
SourceFiles
    tetPolyMeshMapperFaceDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyMeshMapperFaceDecomp_H
#define tetPolyMeshMapperFaceDecomp_H

#include "tetPointMapperFaceDecomp.H"
#include "pointMapper.H"
#include "faceMapper.H"
#include "cellMapper.H"
#include "tetPolyBoundaryMapperFaceDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMeshFaceDecomp;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                 Class tetPolyMeshMapperFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyMeshMapperFaceDecomp
{
    // Private data

        //- Reference to mesh
        const tetPolyMeshFaceDecomp& mesh_;

        //- Point mapper
        pointMapper pointMap_;

        //- Face mapper
        faceMapper faceMap_;

        //- Cell mapper
        cellMapper cellMap_;

        //- Point mapper (needs to be re-done to insert faces and cells)
        tetPointMapperFaceDecomp tetPointMap_;

        //- Boundary mapper
        tetPolyBoundaryMapperFaceDecomp boundaryMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyMeshMapperFaceDecomp(const tetPolyMeshMapperFaceDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyMeshMapperFaceDecomp&);


public:

    // Constructors

        //- Construct from components
        tetPolyMeshMapperFaceDecomp
        (
            const tetPolyMeshFaceDecomp& mesh,
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
        const tetPolyMeshFaceDecomp& mesh() const
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
        const tetPolyBoundaryMapperFaceDecomp& boundaryMap() const
        {
            return boundaryMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
