/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    polyTopoChange

Description
    Class accumulates information on how to perform mesh refinement.
    Once the refinement request is completed, polyMesh modifies the mesh
    topology.

SourceFiles
    polyTopoChange.C

\*---------------------------------------------------------------------------*/

#ifndef polyTopoChange_H
#define polyTopoChange_H

#include "DynamicList.H"
#include "HashSet.H"
#include "polyAddPoint.H"
#include "polyModifyPoint.H"
#include "polyRemovePoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyRemoveFace.H"
#include "polyAddCell.H"
#include "polyModifyCell.H"
#include "polyRemoveCell.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;

/*---------------------------------------------------------------------------*\
                       Class polyTopoChange Declaration
\*---------------------------------------------------------------------------*/

class polyTopoChange
{
    // Private data

        //- Reference to mesh to be refined
        const polyMesh& mesh_;

        //- Points to add
        DynamicList<polyAddPoint> addedPoints_;

        //- Points to modify
        DynamicList<polyModifyPoint> modifiedPoints_;

        //- Points to remove
        labelHashSet removedPoints_;

        //- Faces to add
        DynamicList<polyAddFace> addedFaces_;

        //- Faces to modify
        DynamicList<polyModifyFace> modifiedFaces_;

        //- Faces to remove
        labelHashSet removedFaces_;

        //- Number of cells to add
        DynamicList<polyAddCell> addedCells_;

        //- Cells to modify
        DynamicList<polyModifyCell> modifiedCells_;

        //- Cells to remove
        labelHashSet removedCells_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        polyTopoChange(const polyTopoChange&);

        //- Disallow default bitwise assignment
        void operator=(const polyTopoChange&);


    // Private static data

        //- Minimum dynamic list size for object insertion
        static const label minListSize;

        //- Estimated fraction of removed points
        static const label pointFraction;

        //- Estimated fraction of removed faces
        static const label faceFraction;

        //- Estimated fraction of removed cells
        static const label cellFraction;


public:

    static debug::debugSwitch debug;


    // Constructors

        //- Construct from mesh reference
        polyTopoChange(const polyMesh&);


    // Destructor

        ~polyTopoChange();


    // Member Functions

        // Definition of topological change

            //- Set topological action
            label setAction(const topoAction& action);


        // Topology morph data

            //- Is point removed?
            inline bool pointRemoved(const label pointI) const;

            //- Is face removed?
            inline bool faceRemoved(const label faceI) const;

            //- Is cell removed?
            inline bool cellRemoved(const label cellI) const;

            //- Point balance (added - removed)
            label pointBalance() const
            {
                return
                    addedPoints_.size()
                  + modifiedPoints_.size()
                  - removedPoints_.size();
            }

            //- Face balance (added - removed)
            label faceBalance() const
            {
                return
                    addedFaces_.size()
                  + modifiedFaces_.size()
                  - removedFaces_.size();
            }

            //- Cell balance (added - removed)
            label cellBalance() const
            {
                return addedCells_.size() - removedCells_.size();
            }

            //- Added points
            const DynamicList<polyAddPoint>& addedPoints() const
            {
                return addedPoints_;
            }

            //- Modified points
            const DynamicList<polyModifyPoint>& modifiedPoints() const
            {
                return modifiedPoints_;
            }

            //- Map of removed points
            const labelHashSet& removedPoints() const
            {
                return removedPoints_;
            }

            //- Added faces
            const DynamicList<polyAddFace>& addedFaces() const
            {
                return addedFaces_;
            }

            //- Modified faces
            const DynamicList<polyModifyFace>& modifiedFaces() const
            {
                return modifiedFaces_;
            }

            //- Map of removed faces
            const labelHashSet& removedFaces() const
            {
                return removedFaces_;
            }

            //- Added cells
            const DynamicList<polyAddCell>& addedCells() const
            {
                return addedCells_;
            }

            //- Modified cells
            const DynamicList<polyModifyCell>& modifiedCells() const
            {
                return modifiedCells_;
            }

            //- Map of removed cells
            const labelHashSet& removedCells() const
            {
                return removedCells_;
            }

            //- Check for consistency and report status
            //  Returns false for no error.
            bool check() const;

            //- Clear all contents
            void clear();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "polyTopoChangeI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
