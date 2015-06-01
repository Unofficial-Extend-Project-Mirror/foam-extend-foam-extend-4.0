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
    tetMeshGenerator

Description
    Creates cartesian mesh from the octree

SourceFiles
    tetMeshGenerator.C

\*---------------------------------------------------------------------------*/

#ifndef tetMeshGenerator_H
#define tetMeshGenerator_H

#include "polyMeshGen.H"
#include "IOdictionary.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurf;
class meshOctree;
class Time;

/*---------------------------------------------------------------------------*\
                           Class tetMeshGenerator Declaration
\*---------------------------------------------------------------------------*/

class tetMeshGenerator
{
    // Private data
        //- reference to Time
        const Time& runTime_;

        //- pointer to the surface
        const triSurf* surfacePtr_;

        //- IOdictionary containing information about cell sizes, etc..
        IOdictionary meshDict_;

        //- pointer to the octree
        meshOctree* octreePtr_;

        //- mesh
        polyMeshGen mesh_;

    // Private member functions
        //- create cartesian mesh
        void createTetMesh();

        //- prepare mesh surface
        void surfacePreparation();

        //- map mesh to the surface and untangle surface
        void mapMeshToSurface();

        //- capture edges and corners
        void mapEdgesAndCorners();

        //- optimise surface mesh
        void optimiseMeshSurface();

        //- add boundary layers
        void generateBoudaryLayers();

        //- mesh optimisation
        void optimiseFinalMesh();

        //- refine boundary layers
        void refBoundaryLayers();

        //- replace boundaries
        void replaceBoundaries();

        //- renumber the mesh
        void renumberMesh();

        //- generate mesh
        void generateMesh();

        //- Disallow default bitwise copy construct
        tetMeshGenerator(const tetMeshGenerator&);

        //- Disallow default bitwise assignment
        void operator=(const tetMeshGenerator&);

public:

    // Constructors

        //- Construct from time
        tetMeshGenerator(const Time&);

    // Destructor

        ~tetMeshGenerator();


    // Member Functions

        //- write the mesh
        void writeMesh() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
