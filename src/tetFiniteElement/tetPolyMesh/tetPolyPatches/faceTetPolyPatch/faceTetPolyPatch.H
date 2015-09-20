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
    faceTetPolyPatch

Description
        A tetPolyPatch based on a polyPatch

SourceFiles
    faceTetPolyPatch.C
    faceTetPolyPatchM.C
    newTetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef faceTetPolyPatch_H
#define faceTetPolyPatch_H

#include "tetPolyPatch.H"
#include "polyPatch.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
                  Class faceTetPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class faceTetPolyPatch
:
    public tetPolyPatch
{
    // Private data

        //- Index of polyPatch in boundary
        //  Needs to be stored as faceTetPolyPatch does not hold a reference to
        //  a polyPatch in order to enable topological changes
        const label boundaryIndex_;

        //- Number of points
        label size_;


    // Demand driven private data

        //- Mesh points
        mutable labelList* meshPointsPtr_;

        //- Local points
        mutable pointField* localPointsPtr_;

        //- Point unit normals
        mutable vectorField* pointNormalsPtr_;

        //- Local edges indices
        mutable labelList* localEdgeIndicesPtr_;

        //- Cut edges indices, i.e. the ones originating from the
        //  points on the patch but connected to points which are
        //  internal (not on the patch)
        mutable labelList* cutEdgeIndicesPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceTetPolyPatch(const faceTetPolyPatch&);

        //- Disallow default bitwise assignment
        void operator=(const faceTetPolyPatch&);

        //- Clear out
        void clearOut();


protected:

    // Member functions to calculate demand driven data

        //- Calculate mesh points
        labelList calcMeshPoints(const primitivePatch& p) const;

        //- Calculate local points
        void calcLocalPoints() const;

        //- Calculate point normals
        void calcPointNormals() const;

        //- Calculate cut edge indices
        labelList calcCutEdgeIndices(const primitivePatch& p) const;

        //- Calculate local edges indices given a primitive patch
        labelList calcLocalEdgesIndices(const primitivePatch& p) const;


public:

    typedef tetPolyBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceTetPolyPatch,
            polyPatch,
            (const polyPatch& patch, const tetPolyBoundaryMesh& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch
        faceTetPolyPatch
        (
            const polyPatch&,
            const tetPolyBoundaryMesh&
        );


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<faceTetPolyPatch> New
        (
            const polyPatch&,
            const tetPolyBoundaryMesh&
        );


    // Destructor

        virtual ~faceTetPolyPatch();


    // Member Functions

        //- Return the polyPatch
        const polyPatch& patch() const;

        //- Return name
        virtual const word& name() const
        {
            return patch().name();
        }

        //- Return size
        virtual label size() const
        {
            return size_;
        }

        //- Return number of faces
        virtual label nFaces() const
        {
            return patch().size();
        }

        //- Return the index of this patch in the tetPolyBoundaryMesh
        virtual label index() const
        {
            return boundaryIndex_;
        }


        // Access functions for demand driven data

            //- Return mesh points
            virtual const labelList& meshPoints() const;

            //- Return pointField of points in patch
            virtual const pointField& localPoints() const;

            //- Return point unit normals
            virtual const vectorField& pointNormals() const;

            //- Face decomposition into triangle
            //  Used for implementation of FEM boundary conditions
            virtual triFaceList faceTriangles(const label faceID) const;

            //- Return complete list of faces. All are triangles
            virtual faceList triFaces() const;

            //- Return list of edge indices for edges local to the patch
            //  (i.e. connecting points within the patch)
            virtual const labelList& localEdgeIndices() const;

            //- Return list of edge indices for cut edges
            //  (i.e. connecting points within the patch ot points outside it)
            virtual const labelList& cutEdgeIndices() const;


        // Update topology
        virtual void updateMesh();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
