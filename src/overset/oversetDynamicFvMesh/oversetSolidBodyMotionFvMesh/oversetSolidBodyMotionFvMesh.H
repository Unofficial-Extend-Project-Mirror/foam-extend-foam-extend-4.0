/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
