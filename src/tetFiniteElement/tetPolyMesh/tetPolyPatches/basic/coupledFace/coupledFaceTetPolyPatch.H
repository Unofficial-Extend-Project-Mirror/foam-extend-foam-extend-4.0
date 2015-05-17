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
    coupledFaceTetPolyPatch

Description
    Coupled patch.  Used as the base class for processor
    and cyclic tetPolyPatches

SourceFiles
    coupledFaceTetPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef coupledFaceTetPolyPatch_H
#define coupledFaceTetPolyPatch_H

#include "lduInterface.H"
#include "faceTetPolyPatch.H"
#include "coupledPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class tetPolyBoundaryMesh;

/*---------------------------------------------------------------------------*\
               Class coupledFaceTetPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class coupledFaceTetPolyPatch
:
    public lduInterface,
    public faceTetPolyPatch
{
    // Private data

        const coupledPolyPatch& coupledPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        coupledFaceTetPolyPatch
        (
            const coupledFaceTetPolyPatch&
        );

        //- Disallow default assignment
        void operator=(const coupledFaceTetPolyPatch&);


protected:

    // Demand driven private data

        //- The set of labels of the tetPolyPatch points which are
        //  non-global, i.e. present in this coupled
        mutable labelList* nonGlobalPatchPointsPtr_;

        //- Mesh points
        mutable labelList* meshPointsPtr_;


        // Construction of demand-driven data

            //- Calculate mesh points
            virtual void calcMeshPoints() const = 0;


public:

    typedef tetPolyBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(coupledPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        coupledFaceTetPolyPatch
        (
            const polyPatch& patch,
            const tetPolyBoundaryMesh& bm
        );


    // Destructor

        virtual ~coupledFaceTetPolyPatch();


    // Member Functions

        // Access

            //- Return size
            virtual label size() const
            {
                return meshPoints().size();
            }

            //- Return true because this patch is coupled
            virtual bool coupled() const
            {
                return true;
            }


        // Access functions for demand driven data

            //- Return the set of labels of the tetPolyPatch points which are
            //  non-global, i.e. present in this coupled
            virtual const labelList& nonGlobalPatchPoints() const;

            //- Return mesh points
            virtual const labelList& meshPoints() const;

            //- Return pointField of points in patch.  Not impelemented.
            virtual const pointField& localPoints() const;

            //- Return point unit normals.  Not impelemented.
            virtual const vectorField& pointNormals() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
