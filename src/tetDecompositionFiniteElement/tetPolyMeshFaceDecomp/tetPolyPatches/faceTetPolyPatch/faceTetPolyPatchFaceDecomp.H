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
    faceTetPolyPatchFaceDecomp

Description
        A tetPolyPatch based on a polyPatch

SourceFiles
    faceTetPolyPatchFaceDecomp.C
    faceTetPolyPatchFaceDecompM.C
    newTetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef faceTetPolyPatchFaceDecomp_H
#define faceTetPolyPatchFaceDecomp_H

#include "tetPolyPatchFaceDecomp.H"
#include "polyPatch.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class tetPolyBoundaryMeshFaceDecomp;

/*---------------------------------------------------------------------------*\
                  Class faceTetPolyPatchFaceDecomp Declaration
\*---------------------------------------------------------------------------*/

class faceTetPolyPatchFaceDecomp
:
    public tetPolyPatchFaceDecomp
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
        faceTetPolyPatchFaceDecomp(const faceTetPolyPatchFaceDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const faceTetPolyPatchFaceDecomp&);

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

    typedef tetPolyBoundaryMeshFaceDecomp BoundaryMesh;


    //- Runtime type information
    TypeName("patch");

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceTetPolyPatchFaceDecomp,
            polyPatch,
            (const polyPatch& patch, const tetPolyBoundaryMeshFaceDecomp& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch
        faceTetPolyPatchFaceDecomp
        (
            const polyPatch&,
            const tetPolyBoundaryMeshFaceDecomp&
        );


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<faceTetPolyPatchFaceDecomp> New
        (
            const polyPatch&,
            const tetPolyBoundaryMeshFaceDecomp&
        );


    // Destructor

        virtual ~faceTetPolyPatchFaceDecomp();


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

        //- Return the index of this patch in the tetPolyBoundaryMeshFaceDecomp
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
