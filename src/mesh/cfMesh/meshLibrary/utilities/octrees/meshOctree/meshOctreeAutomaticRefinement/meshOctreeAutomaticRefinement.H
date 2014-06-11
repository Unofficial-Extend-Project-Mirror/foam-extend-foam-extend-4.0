/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Class
    meshOctreeAutomaticRefinement

Description
    Automatic octree for mesh generation. This class checks if the current
    octree resolution is enough to resolve all surface features

SourceFiles
    meshOctreeAutomaticRefinement.C

\*---------------------------------------------------------------------------*/

#ifndef meshOctreeAutomaticRefinement_H
#define meshOctreeAutomaticRefinement_H

#include "meshOctreeModifier.H"
//#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations
class triSurfaceCurvatureEstimator;
class triSurfacePartitioner;
class meshOctreeAddressing;

/*---------------------------------------------------------------------------*\
                Class meshOctreeAutomaticRefinement Declaration
\*---------------------------------------------------------------------------*/

class meshOctreeAutomaticRefinement
{
    // Private data
        //- reference to meshOctree
        meshOctree& octree_;

        //- const reference to mesh dictionary
        const IOdictionary& meshDict_;

        //- use DATA boxes
        bool useDATABoxes_;
    
        //- a flag for activating paired refinement
        bool hexRefinement_;

    // Private member functions
        //- meshOctreeAddressing pointer
        mutable meshOctreeAddressing* octreeAddressingPtr_;

        //- calculate octree addresing
        void createOctreeAddressing() const;
        const meshOctreeAddressing& octreeAddressing() const;

        //- surface partitioner pointer
        mutable triSurfacePartitioner* partitionerPtr_;

        //- calculate partitionerPtr_
        void createSurfacePartitioner() const;
        const triSurfacePartitioner& partitioner() const;

        //- curvature estimator pointer
        mutable triSurfaceCurvatureEstimator* curvaturePtr_;

           //- calculate curvaturePtr_
        void createCurvatureEstimator() const;
        const triSurfaceCurvatureEstimator& curvature() const;

    // Functions and data performing automatic refinement
        //- maximum allowed refinement level
        direction maxRefLevel_;

        //- set maxRefLevel_
        void setMaxRefLevel();

        //- refine selected boxes
        void refineSelectedBoxes(List<direction>&, labelLongList&);

        //- refine DATA boxes based on curvature
        bool curvatureRefinement();
        bool refineBasedOnCurvature(List<direction>&, const labelLongList&);

        //- refine DATA boxes which contain two or more independent parts
        //- of the surface in their vicinity
        bool proximityRefinement();
        bool refineBasedOnProximityTests
        (
            List<direction>&,
            const labelLongList&
        );

        //- refine DATA boxes containing 2 or more corners, 2 or more distinct
        //- edges, etc.
        bool refineBasedOnContainedPartitions
        (
            List<direction>&,
            const labelLongList&
        );

        //- refine boxes based on the number of contained surface corners
        bool refineBasedOnContainedCorners
        (
            List<direction>&,
            const labelLongList&
        );

    // Private copy constructor
        //- Disallow default bitwise copy construct
        meshOctreeAutomaticRefinement(const meshOctreeAutomaticRefinement&);

        //- Disallow default bitwise assignment
        void operator=(const meshOctreeAutomaticRefinement&);

public:

    // Constructors

        //- Construct from meshOctree, dictionary and boolean
        meshOctreeAutomaticRefinement
        (
            meshOctree& mo,
            const IOdictionary& dict,
            bool useDATABoxes = false
        );

    // Destructor

        ~meshOctreeAutomaticRefinement();
        
    // Member Functions
        //- activate hex refinement
        void activateHexRefinement();

        //- refine octree boxes which are too coarse and cannot capture
        //- surface features
        void automaticRefinement();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
