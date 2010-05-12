// The FOAM Project // File: layerAR.H
/*
-------------------------------------------------------------------------------
 =========         | Class Interface
 \\      /         |
  \\    /          | Name:   layerAR
   \\  /           | Family: engine
    \\/            |
    F ield         | FOAM version: 2.3
    O peration     |
    A and          | Copyright (C) 1991-2004 Nabla Ltd.
    M anipulation  |          All Rights Reserved.
-------------------------------------------------------------------------------
CLASS
    layerAR

DESCRIPTION
    Engine mesh class with dynamic layering on the piston patch. 

C++
    listOfSourceFiles
        addlayerARModifiers.C
        layerAR.C
    endListOfSourceFiles

AUTHOR
	
*/
// ------------------------------------------------------------------------- //

#ifndef layerAR_H
#define layerAR_H

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
                           Class layerAR Declaration
\*---------------------------------------------------------------------------*/

class layerAR
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        simpleEnginePiston piston_;
        
        //- Layering-to-deformation switch in crank-angle degrees
        scalar deformSwitch_;
        
        //- Tolerance used when the piston faceZone is created
        scalar delta_;

        //- Distance from the piston patch
        scalar offSet_;

        //- Piston Position
        scalar pistonPosition_;

        //- Virtual piston position (pistonPosition + offSet)
        scalar virtualPistonPosition_;

        //- deckHeight
        scalar deckHeight_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        layerAR(const layerAR&);

        //- Disallow default bitwise assignment
        void operator=(const layerAR&);

        //- Make layering modifiers live
        void makeLayersLive();

        //- Check if all patches exist, then calculate virtualPistonPosition,
        //- pistonPosition and deckHeight for the first time
        void checkAndCalculate();

        //- Calculate the virtualPistonPosition,
        void setVirtualPistonPosition();

        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers();

public:

    //- Runtime type information
    TypeName("layerAR");


    // Constructors

        //- Construct from database
        explicit layerAR(const IOobject& io);


    // Destructor

        virtual ~layerAR();


    // Member Functions

        const scalar& delta() const
        {
            return delta_;
        }

        scalar& delta()
        {
            return delta_;
        }

        const scalar& offSet() const
        {
            return offSet_;
        }

        scalar& offSet()
        {
            return offSet_;
        }

        const scalar& pistonPosition() const
        {
            return pistonPosition_;
        }

        scalar& pistonPosition()
        {
            return pistonPosition_;
        }

        const scalar& virtualPistonPosition() const
        {
            return virtualPistonPosition_;
        }

        scalar& virtualPistonPosition()
        {
            return virtualPistonPosition_;
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

        //- Return true for mesh deformation mode
        bool deformation() const
        {
            return
                engTime().thetaRevolution() > -deformSwitch_
             && engTime().thetaRevolution() < deformSwitch_;
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
