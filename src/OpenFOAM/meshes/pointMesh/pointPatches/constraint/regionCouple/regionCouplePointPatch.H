/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    regionCouplePointPatch

Description
    Dummy regionCouple patch for post-processing.  No functionality built in!

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    regionCouplePointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef regionCouplePointPatch_H
#define regionCouplePointPatch_H

#include "coupledFacePointPatch.H"
#include "regionCouplePolyPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class regionCouplePointPatch Declaration
\*---------------------------------------------------------------------------*/

class regionCouplePointPatch
:
    public coupledFacePointPatch
{
    // Private data

        //- Local reference cast into the cyclic patch
        const regionCouplePolyPatch& regionCouplePolyPatch_;


    // Private Member Functions

        //- Disallow default construct as copy
        regionCouplePointPatch(const regionCouplePointPatch&);

        //- Disallow default assignment
        void operator=(const regionCouplePointPatch&);


protected:

    // Protected Member Functions

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
    TypeName(regionCouplePolyPatch::typeName_());


    // Constructors

        //- Construct from components
        regionCouplePointPatch
        (
            const polyPatch& patch,
            const pointBoundaryMesh& bm
        );


    // Destructor

        virtual ~regionCouplePointPatch();


    // Member Functions

        // Access

            //- Return true if coupled
            virtual bool coupled() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
