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

        //- Immersed boundary motion control function
        PtrList<movingImmersedBoundary> ibMotions_;

        //- Current mesh motion time index
        label curTimeIndex_;


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
