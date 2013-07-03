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
    mixingPlanePointPatch

Description
    Dummy mixingPlane patch for post-processing.  No functionality built in!

Author
    Martin Beaudoin, Hydro-Quebec, 2009.  All rights reserved

Contributor
    Hrvoje Jasak, Wikki Ltd.

SourceFiles
    mixingPlanePointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef mixingPlanePointPatch_H
#define mixingPlanePointPatch_H

#include "coupledFacePointPatch.H"
#include "mixingPlanePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class mixingPlanePointPatch Declaration
\*---------------------------------------------------------------------------*/

class mixingPlanePointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const mixingPlanePolyPatch& mixingPlanePolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        mixingPlanePointPatch(const mixingPlanePointPatch&);

        //- Disallow default assignment
        void operator=(const mixingPlanePointPatch&);


        //- Initialise the calculation of the patch geometry
        virtual void initGeometry()
        {}

        //- Calculate mesh points
        virtual void calcGeometry()
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
    TypeName(mixingPlanePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        mixingPlanePointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~mixingPlanePointPatch();


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
