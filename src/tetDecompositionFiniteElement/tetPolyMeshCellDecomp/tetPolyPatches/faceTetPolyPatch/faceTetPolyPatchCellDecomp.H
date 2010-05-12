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
    faceTetPolyPatchCellDecomp

Description
    A tetPolyPatch based on a polyPatch

SourceFiles
    faceTetPolyPatchCellDecomp.C
    faceTetPolyPatchCellDecompM.C
    newTetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef faceTetPolyPatchCellDecomp_H
#define faceTetPolyPatchCellDecomp_H

#include "tetPolyPatchCellDecomp.H"
#include "polyPatch.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class faceTetPolyPatchCellDecomp Declaration
\*---------------------------------------------------------------------------*/

class faceTetPolyPatchCellDecomp
:
    public tetPolyPatchCellDecomp
{
    // Private data

        //- Index of polyPatch in boundary
        //  Needs to be stored as faceTetPolyPatch does not hold a reference to
        //  a polyPatch in order to enable topological changes
        const label boundaryIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        faceTetPolyPatchCellDecomp(const faceTetPolyPatchCellDecomp&);

        //- Disallow default bitwise assignment
        void operator=(const faceTetPolyPatchCellDecomp&);


public:

    typedef tetPolyBoundaryMeshCellDecomp BoundaryMesh;


    //- Runtime type information
    TypeName(polyPatch::typeName_());

    // Declare run-time constructor selection tables

        declareRunTimeSelectionTable
        (
            autoPtr,
            faceTetPolyPatchCellDecomp,
            polyPatch,
            (const polyPatch& patch, const tetPolyBoundaryMeshCellDecomp& bm),
            (patch, bm)
        );


    // Constructors

        //- Construct from polyPatch
        faceTetPolyPatchCellDecomp
        (
            const polyPatch&,
            const tetPolyBoundaryMeshCellDecomp&
        );


    // Selectors

        //- Return a pointer to a new patch created on freestore from polyPatch
        static autoPtr<faceTetPolyPatchCellDecomp> New
        (
            const polyPatch&,
            const tetPolyBoundaryMeshCellDecomp&
        );


    // Destructor

        virtual ~faceTetPolyPatchCellDecomp()
        {}


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
            return patch().nPoints();
        }

        //- Return number of faces
        virtual label nFaces() const
        {
            return patch().size();
        }

        //- Return the index of this patch in the tetPolyBoundaryMeshCellDecomp
        virtual label index() const
        {
            return boundaryIndex_;
        }


        // Access functions for demand driven data

            //- Return mesh points
            virtual const labelList& meshPoints() const
            {
                return patch().meshPoints();
            }

            //- Return pointField of points in patch
            virtual const pointField& localPoints() const
            {
                return patch().localPoints();
            }

            //- Return point unit normals
            virtual const vectorField& pointNormals() const
            {
                return patch().pointNormals();
            }

            //- Face decomposition into triangle
            //  Used for implementation of FEM boundary conditions
            virtual triFaceList faceTriangles(const label faceID) const;

            //- Return complete list of faces. All are triangles
            virtual faceList triFaces() const;


        // Update topology - not required
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
