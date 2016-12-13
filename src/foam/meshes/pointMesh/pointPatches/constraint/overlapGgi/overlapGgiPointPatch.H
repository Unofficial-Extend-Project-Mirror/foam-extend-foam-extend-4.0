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
    overlapGgiPointPatch

Description
    Dummy overlapGgi patch for post-processing.  No functionality built in!

Author:
    Hrvoje Jasak, Wikki Ltd.  All rigths reserved.

SourceFiles
    overlapGgiPointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiPointPatch_H
#define overlapGgiPointPatch_H

#include "coupledFacePointPatch.H"
#include "overlapGgiPolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class overlapGgiPointPatch Declaration
\*---------------------------------------------------------------------------*/

class overlapGgiPointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the overlap patch
        const overlapGgiPolyPatch& overlapGgiPolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        overlapGgiPointPatch(const overlapGgiPointPatch&);

        //- Disallow default assignment
        void operator=(const overlapGgiPointPatch&);


protected:

    // Protected Member Functions

        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate mesh points
        virtual void calcGeometry()
        {}

        //- Initialise the patches for moving points
        virtual void initMovePoints(const pointField&)
        {}

        //- Correct patches after moving points
        virtual void movePoints()
        {}

        //- Initialise the update of the patch topology
        virtual void initUpdateMesh();

        //- Update of the patch topology
        virtual void updateMesh();


public:

    typedef pointBoundaryMesh BoundaryMesh;


    //- Runtime type information
    TypeName(overlapGgiPolyPatch::typeName_());


    // Constructors

        //- Construct from components
        overlapGgiPointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~overlapGgiPointPatch();


    // Member Functions

        // Access

            //- Return true because this patch is coupled
            virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
