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
    Foam::immersedBoundaryPolyPatch

Description
    Immersed boundary patch

Author
    Zeljko Tukovic
    Reorganisation by Hrvoje Jasak

SourceFiles
    immersedBoundaryPolyPatch.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundaryPolyPatch_H
#define immersedBoundaryPolyPatch_H

#include "polyPatch.H"
#include "triSurfaceMesh.H"
#include "triSurfaceTools.H"
#include "triSurfaceSearch.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class immersedBoundaryPolyPatch Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundaryPolyPatch
:
    public polyPatch
{
    // Private data

        //- Triangular surface representing immersed boundary.
        //  Name of tri surface will be identical to the name of the patch
        triSurfaceMesh ibMesh_;

        //- Internal or external flow calculation
        Switch internalFlow_;


    // Demand-driven data

        //- Triangular surface search algorithm
        mutable triSurfaceSearch* triSurfSearchPtr_;

        //- Moving immersed boundary
        bool movingIb_;


    // Private Member Functions

        // Storage management

            //- Clear all demand-driven data
            void clearOut();


        // Make demand-driven data

            //- Make triSurface search algorithm
            void makeTriSurfSearch() const;


protected:

    // Protected Member Functions

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints(const pointField& p);


public:

    //- Runtime type information
    TypeName("immersedBoundary");


    // Constructors

        //- Construct from components
        immersedBoundaryPolyPatch
        (
            const word& name,
            const label size,
            const label start,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct from dictionary
        immersedBoundaryPolyPatch
        (
            const word& name,
            const dictionary& dict,
            const label index,
            const polyBoundaryMesh& bm
        );

        //- Construct as copy, resetting the boundary mesh
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch&,
            const polyBoundaryMesh&
        );

        //- Construct given the original patch and resetting the
        //  face list and boundary mesh information
        immersedBoundaryPolyPatch
        (
            const immersedBoundaryPolyPatch& pp,
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        );

        //- Construct and return a clone, resetting the boundary mesh
        virtual autoPtr<polyPatch> clone(const polyBoundaryMesh& bm) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch(*this, bm)
            );
        }

        //- Construct and return a clone, resetting the face list
        //  and boundary mesh
        virtual autoPtr<polyPatch> clone
        (
            const polyBoundaryMesh& bm,
            const label index,
            const label newSize,
            const label newStart
        ) const
        {
            return autoPtr<polyPatch>
            (
                new immersedBoundaryPolyPatch
                (
                    *this,
                    bm,
                    index,
                    newSize,
                    newStart
                )
            );
        }


    //- Destructor
    virtual ~immersedBoundaryPolyPatch();


    // Member Functions

        // Access

            //- Return immersed boundary surface mesh
            const triSurfaceMesh& ibMesh() const
            {
                return ibMesh_;
            }

            //- Return true if solving for flow inside the immersed boundary
            bool internalFlow() const
            {
                return internalFlow_;
            }

            //- Return triSurface search object
            const triSurfaceSearch& triSurfSearch() const;

            //- Return true if immersed boundary is moving
            bool movingIb() const
            {
                return movingIb_;
            }


        // Edit

            //- Correct patches after moving points
            void moveTriSurfacePoints(const pointField& p);



        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
