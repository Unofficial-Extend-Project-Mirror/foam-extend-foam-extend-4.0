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
    tetPolyMeshMapperCellDecomp

Description
    Class holds all the necessary information for mapping fields associated
    with tetPolyMeshCellDecomp.
SourceFiles
    tetPolyMeshMapperCellDecomp.C

\*---------------------------------------------------------------------------*/

#ifndef tetPolyMeshMapperCellDecomp_H
#define tetPolyMeshMapperCellDecomp_H

#include "tetPointMapperCellDecomp.H"
#include "pointMapper.H"
#include "cellMapper.H"
#include "tetPolyBoundaryMapperCellDecomp.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyMeshCellDecomp;
class mapPolyMesh;

/*---------------------------------------------------------------------------*\
                 Class tetPolyMeshMapperCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class tetPolyMeshMapperCellDecomp
{
    // Private data

        //- Reference to mesh
        const tetPolyMeshCellDecomp& mesh_;

        //- Point mapper
        pointMapper pointMap_;

        //- Cell mapper
        cellMapper cellMap_;

        //- Point mapper (needs to be re-done to insert cells)
        tetPointMapperCellDecomp tetPointMap_;

        //- Boundary mapper
        tetPolyBoundaryMapperCellDecomp boundaryMap_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        tetPolyMeshMapperCellDecomp(const tetPolyMeshMapperCellDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const tetPolyMeshMapperCellDecomp&);


public:

    // Constructors

        //- Construct from components
        tetPolyMeshMapperCellDecomp
        (
            const tetPolyMeshCellDecomp& mesh,
            const mapPolyMesh& meshMap
        )
        :
            mesh_(mesh),
            pointMap_(meshMap),
            cellMap_(meshMap),
            tetPointMap_(mesh, meshMap, pointMap_, cellMap_),
            boundaryMap_(mesh, meshMap)
        {}


    // Member Functions

        //- Return mesh
        const tetPolyMeshCellDecomp& mesh() const
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
        const tetPolyBoundaryMapperCellDecomp& boundaryMap() const
        {
            return boundaryMap_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
