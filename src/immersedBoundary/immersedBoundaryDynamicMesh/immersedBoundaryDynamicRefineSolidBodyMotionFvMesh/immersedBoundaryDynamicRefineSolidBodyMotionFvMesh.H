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
    Foam::immersedBoundaryDynamicRefineSolidBodyMotionFvMesh

Description
    A fvMesh with built-in refinement of arbitrary polyhedral cells.

    Determines which cells to refine/unrefine and does all in update().

SourceFiles
    immersedBoundaryDynamicRefineSolidBodyMotionFvMesh.C

Author
    Vuko Vukcevic, Wikki Ltd.  All rights reserved
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

Notes
    Dynamcially refining immersed boundary mesh with solid body motion
    of the immersed boundary

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryDynamicRefineSolidBodyMotionFvMesh_H
#define immersedBoundaryDynamicRefineSolidBodyMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "polyRef.H"
#include "PackedBoolList.H"
#include "movingImmersedBoundary.H"
#include "triSurfaceDistance.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
     Class immersedBoundaryDynamicRefineSolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryDynamicRefineSolidBodyMotionFvMesh
:
    public dynamicFvMesh
{

    // Private data

        //- Immersed boundary motion control function
        PtrList<movingImmersedBoundary> ibMotions_;

        //- Mesh cutting engine
        polyRef meshCutter_;

        //- Current mesh motion time index
        label curTimeIndex_;

        //- Number of refinement/unrefinement steps done so far.
        label nRefinementIterations_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        immersedBoundaryDynamicRefineSolidBodyMotionFvMesh
        (
            const immersedBoundaryDynamicRefineSolidBodyMotionFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const immersedBoundaryDynamicRefineSolidBodyMotionFvMesh&
        );


        //- Count set/unset elements in packedlist.
        static label count(const PackedBoolList&, const unsigned int);

        //- Refine cells. Update mesh and fields.
        autoPtr<mapPolyMesh> refine(const labelList&);

        //- Unrefine cells. Gets passed in centre points of cells to combine.
        autoPtr<mapPolyMesh> unrefine(const labelList&);


        // Selection of cells to refine and coarsen

            //- Subset candidate cells for refinement
            labelList selectRefineCells
            (
                const label maxCells,
                const label maxRefLevel,
                const scalar refinementDistance,
                PackedBoolList& candidateCell
            ) const;

            //- Select points that can be unrefined
            labelList selectUnrefinePoints
            (
                const scalar coarsenDistance,
                const PackedBoolList& markedCell
            ) const;

            //- Extend markedCell with cell-face-cell
            void extendMarkedCells(PackedBoolList& markedCell) const;


public:

    //- Runtime type information
    TypeName("immersedBoundaryDynamicRefineSolidBodyMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit immersedBoundaryDynamicRefineSolidBodyMotionFvMesh
        (
            const IOobject& io
        );


    //- Destructor
    virtual ~immersedBoundaryDynamicRefineSolidBodyMotionFvMesh();


    // Member Functions

        //- Direct access to the refinement engine
        const polyRef& meshCutter() const
        {
            return meshCutter_;
        }

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
