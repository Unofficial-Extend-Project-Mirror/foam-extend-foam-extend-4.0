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
    Foam::processorPointPatch

Description
    Processor patch boundary needs to be such that the ordering of
    points in the patch is the same on both sides.

    Looking at the creation of the faces on both sides of the processor
    patch they need to be identical on both sides with the normals pointing
    in opposite directions.  This is achieved by calling the reverseFace
    function in the decomposition.  It is therefore possible to re-create
    the ordering of patch points on the slave side by reversing all the
    patch faces of the owner.

SourceFiles
    processorPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef processorPointPatch_H
#define processorPointPatch_H

#include "coupledFacePointPatch.H"
#include "processorPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class processorPointPatch Declaration
\*---------------------------------------------------------------------------*/

class processorPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Reference to processor patch
        const processorPolyPatch& procPolyPatch_;

        //- Mesh points.  On slave they will be reversed
        mutable labelList meshPoints_;


    // Private Member Functions

        //- Disallow default construct as copy
        processorPointPatch(const processorPointPatch&);

        //- Disallow default assignment
        void operator=(const processorPointPatch&);


        //- Initialise the points on this patch which are should also be
        //  on a neighbouring patch but are not part of faces of that patch
        void initPatchPatchPoints();

        //- Calculate the points on this patch which are should also be
        //  on a neighbouring patch but are not part of faces of that patch
        void calcPatchPatchPoints();


protected:

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry();

        //- Calculate the patch geometry
        virtual void calcGeometry();

        //- Correct patches after moving points
        virtual void movePoints()
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    typedef coupledPointPatch CoupledPointPatch;
    typedef pointBoundaryMesh BoundaryMesh;

    //- Runtime type information
    TypeName(processorPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        processorPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~processorPointPatch();


    // Member functions

        //- Return true if running parallel
        virtual bool coupled() const
        {
            if (Pstream::parRun())
            {
                return true;
            }
            else
            {
                return false;
            }
        }

        //- Return processor number
        int myProcNo() const
        {
            return procPolyPatch_.myProcNo();
        }

        //- Return neigbour processor number
        int neighbProcNo() const
        {
            return procPolyPatch_.neighbProcNo();
        }

        //- Is this a master patch
        bool isMaster() const
        {
            return myProcNo() < neighbProcNo();
        }

        //- Is this a slave patch
        bool isSlave() const
        {
            return !isMaster();
        }

        //- Return the underlying processorPolyPatch
        const processorPolyPatch& procPolyPatch() const
        {
            return procPolyPatch_;
        }


        // Access functions for demand driven data

            //- Return mesh points
            virtual const labelList& meshPoints() const;

            //- Return pointField of points in patch.  Not implemented
            virtual const pointField& localPoints() const;

            //- Return point unit normals.  Not implemented
            virtual const vectorField& pointNormals() const;

            //- Return list of edge indices for edges local to the patch
            //  (i.e. connecting points within the patch)
            virtual const labelList& localEdgeIndices() const;

            //- Return list of edge indices for cut edges
            //  (i.e. connecting points within the patch ot points outside it)
            virtual const labelList& cutEdgeIndices() const;


            // Cut edge addressing

                //- Return cut edge owner edge indices
                const labelList& cutEdgeOwnerIndices() const;

                //- Return cut edge owner edge starts
                const labelList& cutEdgeOwnerStart() const;

                //- Return cut edge neighbour edge indices
                const labelList& cutEdgeNeighbourIndices() const;

                //- Return cut edge neighbour edge starts
                const labelList& cutEdgeNeighbourStart() const;

                //- Return doubly cut edge indices
                const labelList& doubleCutEdgeIndices() const;

               //- Return doubly cut edge owner addressing
               //  into current patch
               const labelList& doubleCutOwner() const;

               //- Return doubly cut edge neighbour addressing
               //  into current patch
               const labelList& doubleCutNeighbour() const;

              //- Return cut edge multiplication mask
              const scalarField& ownNeiDoubleMask() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
