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
    Foam::oversetSolidBodyMotionFvMesh

Description
    Solid-body motion of the overset mesh specified by a run-time selectable
    motion function.

SourceFiles
    oversetSolidBodyMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef oversetSolidBodyMotionFvMesh_H
#define oversetSolidBodyMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "movingOversetRegion.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class oversetSolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class oversetSolidBodyMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- Overset region motion control function
        PtrList<movingOversetRegion> motionRegions_;

        //- Reference points which are transformed
        pointIOField undisplacedPoints_;


        // Helper variables

            //- Current time index
            label curTimeIndex_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        oversetSolidBodyMotionFvMesh(const oversetSolidBodyMotionFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const oversetSolidBodyMotionFvMesh&);


public:

    //- Runtime type information
    TypeName("oversetSolidBodyMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit oversetSolidBodyMotionFvMesh(const IOobject& io);


    // Destructor
    virtual ~oversetSolidBodyMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
