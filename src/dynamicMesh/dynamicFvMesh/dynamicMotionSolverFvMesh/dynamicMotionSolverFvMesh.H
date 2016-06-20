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
    Foam::dynamicMotionSolverFvMesh

Description
    Dynamic FV mesh, where a motion solver is used to move the mesh.
    The user specifies motion using the boundary condition on the
    appropriate motion variable.

SourceFiles
    dynamicMotionSolverFvMesh.C

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef dynamicMotionSolverFvMesh_H
#define dynamicMotionSolverFvMesh_H

#include "dynamicFvMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class motionSolver;

/*---------------------------------------------------------------------------*\
                   Class dynamicMotionSolverFvMesh Declaration
\*---------------------------------------------------------------------------*/

class dynamicMotionSolverFvMesh
:
    public dynamicFvMesh
{
    // Private data

        autoPtr<motionSolver> motionPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        dynamicMotionSolverFvMesh(const dynamicMotionSolverFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const dynamicMotionSolverFvMesh&);


public:

    //- Runtime type information
    TypeName("dynamicMotionSolverFvMesh");


    // Constructors

        //- Construct from IOobject
        dynamicMotionSolverFvMesh(const IOobject& io);


    // Destructor

        virtual ~dynamicMotionSolverFvMesh();


    // Member Functions

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
