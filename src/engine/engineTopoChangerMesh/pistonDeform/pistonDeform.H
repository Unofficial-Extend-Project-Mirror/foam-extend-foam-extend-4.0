// The FOAM Project // File: pistonDeform.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   pistonDeform
   \\  /           | Family: engine
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    pistonDeform

DESCRIPTION
    Engine mesh class with dynamic layering on the piston patch. 

C++
    listOfSourceFiles
        addpistonDeformModifiers.C
        pistonDeform.C
    endListOfSourceFiles

AUTHOR
	
*/
// ------------------------------------------------------------------------- //

#ifndef pistonDeform_H
#define pistonDeform_H

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
                           Class pistonDeform Declaration
\*---------------------------------------------------------------------------*/

class pistonDeform
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        simpleEnginePiston piston_;

        //- Piston Position
        scalar pistonPosition_;

        //- deckHeight
        scalar deckHeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        pistonDeform(const pistonDeform&);

        //- Disallow default bitwise assignment
        void operator=(const pistonDeform&);

        //- Check if all patches exist, then calculate virtualPistonPosition,
        //  pistonPosition and deckHeight for the first time
        void checkAndCalculate();

        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers()
        {}

public:

    //- Runtime type information
    TypeName("pistonDeform");


    // Constructors

        //- Construct from database
        explicit pistonDeform(const IOobject& io);


    // Destructor

        virtual ~pistonDeform();

    // Member Functions


        const scalar& pistonPosition() const
        {
            return pistonPosition_;
        }

        scalar& pistonPosition()
        {
            return pistonPosition_;
        }

        const scalar& deckHeight() const
        {
            return deckHeight_;
        }

        scalar& deckHeight()
        {
            return deckHeight_;
        }

        const simpleEnginePiston& piston() const
        {
            return piston_;
        }


        //- Move and morph
        virtual bool update();

        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
