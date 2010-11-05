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
    eMesh

Description
    Mesh class to store edge-based connectivity structures.

SourceFiles
    eMesh.C

\*---------------------------------------------------------------------------*/

#ifndef eMesh_H
#define eMesh_H

#include "Time.H"
#include "polyMesh.H"
#include "edgeIOList.H"
#include "labelIOList.H"
#include "eBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class eMesh Declaration
\*---------------------------------------------------------------------------*/

class eMesh
:
    public objectRegistry
{
    // Private data

        const polyMesh& mesh_;

        //- Ordered edge-list
        mutable edgeIOList edges_;

        //- Boundary mesh
        mutable eBoundaryMesh boundary_;

        // Primitive size data

            //- Number of edges
            mutable label nEdges_;

            //- Number of internal edges
            mutable label nInternalEdges_;

    // Demand-driven data

            //- Mapping between ordered edge-data and regular edges.
            labelList reverseEdgeMap_;

            //- Point-edges
            mutable labelListList* pePtr_;

            //- Edge-points
            mutable labelListList* epPtr_;

            //- Face-edges
            mutable labelListIOList* fePtr_;

            //- Edge-faces
            mutable labelListIOList* efPtr_;

    // Private Member Functions

        //- Helper function to isolate points on triangular faces
        label findIsolatedPoint(const face& f, const edge& e) const;

        //- Helper function to determine the orientation of a triangular face
        label edgeDirection(const face& f, const edge& e) const;

        //- Disallow default bitwise copy construct
        eMesh(const eMesh&);

        //- Disallow default bitwise assignment
        void operator=(const eMesh&);

        // Private member functions to calculate demand driven data

            //- Calculate ordered edges (and edgeFaces)
            void calcOrderedEdgeList();

            //- Calculate ordered edgePoints
            void calcEdgePoints() const;

            //- Calculate point-edges
            void calcPointEdges() const;

            //- Calculate face-edges
            void calcFaceEdges() const;

            //- Calculate edge-faces
            void calcEdgeFaces() const;

            //- Clear geometry
            void clearGeom() const;

            //- Clear addressing
            void clearAddressing() const;

public:

    // Public typedefs

        typedef eMesh Mesh;
        typedef eBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName("eMesh");

    //- Return the mesh sub-directory name (usually "eMesh")
    static word meshSubDir;


    // Constructors

        //- Construct from IOobject and polyMesh reference
        eMesh(const polyMesh& m, const word& subDir = eMesh::meshSubDir);


    // Destructor

        virtual ~eMesh();


    // Member Functions

        // Helpers

            //- Add boundary patches. Constructor helper
            void addEdgePatches(const List<ePatch*> &);


        // Database

            //- Return reference to the mesh database
            virtual const objectRegistry& db() const;

            //- Return the base mesh directory (dbDir())
            fileName meshDir() const;

            //- Return the local mesh directory (dbDir()/meshSubDir)
            fileName meshSubDirectory() const;

            //- Return reference to time
            const Time& time() const;


            //- Mesh size parameters

                label nEdges() const;

                label nInternalEdges() const;

            // Primitive mesh data

                //- Return constant reference to the ordered edge-list
                const edgeList& edges() const;

        // Access

            //- Return constant reference to boundary mesh
            const eBoundaryMesh& boundary() const;


        // Demand-driven data

            //- Return constant reference to the pointEdges list
            const labelListList& pointEdges() const;

            //- Return constant reference to the edgePoints list
            //- Valid for 3D tetrahedral meshes only.
            const labelListList& edgePoints() const;

            //- Return constant reference to the faceEdges list
            const labelListList& faceEdges() const;

            //- Return constant reference to the edgeFaces list
            const labelListList& edgeFaces() const;

        // Reset primitive data
        void resetPrimitives
        (
            edgeList& edges,
            labelListList& faceEdges,
            labelListList& edgeFaces,
            const labelList& patchSizes,
            const labelList& patchStarts,
            const bool reUse,
            const bool storePrimitives
        );

        //- Clear demand-driven data
        void clearOut() const;

        //- Set the instance for mesh files
        void setInstance(const fileName& inst);

        //- Write mesh
        virtual bool write() const;


    // Member Operators

        bool operator!=(const eMesh& m) const;

        bool operator==(const eMesh& m) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
