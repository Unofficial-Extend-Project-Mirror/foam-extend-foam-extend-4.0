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
    Foam::immersedBoundarySolidBodyMotionFvMesh

Description
    Solid-body motion of the immersed boundary mesh specified by a
    run-time selectable motion function for each immersed boundary surface

SourceFiles
    immersedBoundarySolidBodyMotionFvMesh.C

\*---------------------------------------------------------------------------*/

#ifndef immersedBoundarySolidBodyMotionFvMesh_H
#define immersedBoundarySolidBodyMotionFvMesh_H

#include "dynamicFvMesh.H"
#include "dictionary.H"
#include "movingImmersedBoundary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class immersedBoundarySolidBodyMotionFvMesh Declaration
\*---------------------------------------------------------------------------*/

class immersedBoundarySolidBodyMotionFvMesh
:
    public dynamicFvMesh
{
    // Private data

        //- Dictionary of motion control parameters
        dictionary dynamicMeshCoeffs_;

        //- Immersed boundary motion control function
        PtrList<movingImmersedBoundary> ibMotions_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        immersedBoundarySolidBodyMotionFvMesh
        (
            const immersedBoundarySolidBodyMotionFvMesh&
        );

        //- Disallow default bitwise assignment
        void operator=
        (
            const immersedBoundarySolidBodyMotionFvMesh&
        );


public:

    //- Runtime type information
    TypeName("immersedBoundarySolidBodyMotionFvMesh");


    // Constructors

        //- Construct from IOobject
        explicit immersedBoundarySolidBodyMotionFvMesh(const IOobject& io);


    // Destructor
    virtual ~immersedBoundarySolidBodyMotionFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
