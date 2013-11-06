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
    convexSetAlgorithm

Description
    Base class for convexSetAlgorithms

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

SourceFiles
    convexSetAlgorithm.C

\*---------------------------------------------------------------------------*/

#ifndef convexSetAlgorithm_H
#define convexSetAlgorithm_H

#include "Map.H"
#include "label.H"
#include "edgeList.H"
#include "faceList.H"
#include "cellList.H"
#include "boundBox.H"
#include "objectMap.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                    Class convexSetAlgorithm Declaration
\*---------------------------------------------------------------------------*/

class convexSetAlgorithm
{

protected:

    // Protected data

        const label nOldPoints_;

        //- References to old-level connectivity
        //  [Before topo-changes, at old point-positions]
        const polyMesh& mesh_;

        //- References to new-level connectivity
        //  [After topo-changes, at old point-positions]
        const pointField& newPoints_;
        const UList<edge>& newEdges_;
        const UList<face>& newFaces_;
        const UList<cell>& newCells_;
        const UList<label>& newOwner_;
        const UList<label>& newNeighbour_;

        //- Entity parents
        mutable labelList parents_;

        //- Internal data members
        mutable boundBox box_;
        mutable vector refNorm_;
        mutable scalar normFactor_;
        mutable vectorField centres_;
        mutable scalarField weights_;

public:

    //- Constructor

        // Construct from components
        convexSetAlgorithm
        (
            const polyMesh& mesh,
            const pointField& newPoints,
            const UList<edge>& newEdges,
            const UList<face>& newFaces,
            const UList<cell>& newCells,
            const UList<label>& newOwner,
            const UList<label>& newNeighbour
        );

    //- Destructor

        virtual ~convexSetAlgorithm()
        {}

    //- Member functions

        // Dimensions of the algorithm
        virtual label dimension() const = 0;

        // Return true if accumulated weights are consistent
        virtual bool consistent(const scalar tolerance) const;

        // Return the normFactor
        virtual scalar normFactor() const;

        // Normalize stored weights
        virtual void normalize(bool normSum) const;

        // Check whether the bounding box contains the entity
        virtual bool contains(const label index) const = 0;

        // Extract weights and centres to lists
        virtual void populateLists
        (
            labelList& parents,
            vectorField& centres,
            scalarField& weights
        ) const;

        // Compute normFactor
        virtual void computeNormFactor(const label index) const = 0;

        // Compute intersections
        virtual bool computeIntersection
        (
            const label newIndex,
            const label oldIndex,
            const label offset,
            bool output
        ) const = 0;

        // Obtain map weighting factors
        virtual void computeWeights
        (
            const label index,
            const label offset,
            const labelList& mapCandidates,
            const labelListList& oldNeighbourList,
            labelList& parents,
            scalarField& weights,
            vectorField& centres,
            bool output = false
        );

        // Write out connectivity information to disk
        bool write() const;

        // Output an entity as a VTK file
        void writeVTK
        (
            const word& name,
            const label entity,
            const label primitiveType = 3,
            const bool useOldConnectivity = false
        ) const;

        // Output a list of entities as a VTK file
        void writeVTK
        (
            const word& name,
            const labelList& cList,
            const label primitiveType = 3,
            const bool useOldConnectivity = false
        ) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
