/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2006-7 H. Jasak All rights reserved
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
    turboFvMesh

Description
    Simple mixer mesh using an overlapGGI interface

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.

SourceFiles
    turboFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef turboFvMesh_H
#define turboFvMesh_H

#include "dynamicFvMesh.H"
#include "cylindricalCS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class turboFvMesh Declaration
\*---------------------------------------------------------------------------*/

class turboFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Motion dictionary
        dictionary dict_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;

        // - Rotational speed in rotations per minute (rpm)
       mutable scalar rpm_;

        //- Markup field for points.
        mutable vectorField* movingPointsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        turboFvMesh(const turboFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const turboFvMesh&);


        //- Add mixer zones
        void addZonesAndModifiers();

        //- Calculate moving Points
        void calcMovingPoints() const;

        //- Return moving points
        const vectorField& movingPoints() const;


public:

    //- Runtime type information
    TypeName("turboFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit turboFvMesh(const IOobject& io);


    // Destructor

        virtual ~turboFvMesh();


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Update the mesh for both mesh motion
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
