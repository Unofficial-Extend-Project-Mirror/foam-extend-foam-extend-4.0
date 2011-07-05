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
    Simple mixer mesh using an ggi interfaces

    This tool is used to have multiple rotating regions around the same origin
    with different rpms. Creating the cellZones is not implemented
    in this tool.
    The steps to obtain the cellZones and faceZones are:

    1) use regionCellSets utility. With this command you can have different
       cellSets for each region.

    2) run command "setsToZones -noFlipMap".  After this command the
       cellSets are transformed to cellZones.

    4) for each rotating cellZone add an entry in rpm subDict of
       constant/dynamicMeshDict

    5) in parallel you need to create from all rotating coupled-interface
       patches faceZones (if you use a ggi interface between you have already
       created these faceZones), Then you need to specify these faceZones with
       the corresponding rpm in the additional "slider" subDict in
       constant/dynamicMeshDict

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved.
    Fethi Tekin, All rights reserved.
    Oliver Borm, All rights reserved.

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

        //- Markup field for points.
        mutable vectorField* movingPointsPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        turboFvMesh(const turboFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const turboFvMesh&);

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
