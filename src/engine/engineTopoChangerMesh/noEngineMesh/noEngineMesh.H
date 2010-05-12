// The FOAM Project // File: noEngineMesh.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   noEngineMesh
   \\  /           | Family: engine
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    noEngineMesh

DESCRIPTION
    Engine mesh class with dynamic layering on the piston patch. 

C++
    listOfSourceFiles
        addnoEngineMeshModifiers.C
        noEngineMesh.C
    endListOfSourceFiles

AUTHOR
	
*/
// ------------------------------------------------------------------------- //

#ifndef noEngineMesh_H
#define noEngineMesh_H

#include "engineTopoChangerMesh.H"
#include "simpleEnginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class noEngineMesh Declaration
\*---------------------------------------------------------------------------*/

class noEngineMesh
:
    public engineTopoChangerMesh
{
               
    // Private Member Functions

        //- Disallow default bitwise copy construct
        noEngineMesh(const noEngineMesh&);

        //- Disallow default bitwise assignment
        void operator=(const noEngineMesh&);


public:

    //- Runtime type information
    TypeName("noEngineMesh");


    // Constructors

        //- Construct from database
        explicit noEngineMesh(const IOobject& io);


    // Destructor 

        virtual ~noEngineMesh();

    // Member Functions


        //- Move and morph
        virtual bool update()
        {
            return false;
        }
        
        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U){};

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
