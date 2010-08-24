/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2005-2010 Tommaso Lucchini
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
        

Description


SourceFiles

Author
	Tommaso Lucchini
	
\*---------------------------------------------------------------------------*/

#ifndef engineValveSliding_H
#define engineValveSliding_H

#include "engineTopoChangerMesh.H"
#include "enginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "thoboisSlidingValveBank.H"
#include "twoDPointCorrector.H"

#include "faceSet.H"
#include "pointSet.H"
#include "cellSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class verticalValves Declaration
\*---------------------------------------------------------------------------*/

class engineValveSliding
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        enginePiston piston_;
        
        //- vertical valves
        thoboisSlidingValveBank valves_;
                         
        //- Layering-to-deformation switch in crank-angle degrees
        scalar deformSwitch_;
        
        //- Top valve tolerance, when valve lift < tol, valve top layering is disabled
        scalar valveTopTol_;

        //- Piston Position
        scalar pistonPosition_;
        
        //- Virtual position of the bottom of the valves
        scalarField valveBottomPosition_;

        //- Virtual position of the bottom of the valves
        scalarField valvePistonPosition_;

        //- deckHeight
        scalar deckHeight_;
        
        //- min valveHead z
        scalarField minValveZ_;
                        
        //- tolerance for valve poppet faces (used to create a faceZone)
        scalar poppetValveTol_;
        
        //- tolerance for valve poppet faces (used to create a faceZone)
        scalar bottomValveTol_;
        
        //- Motion solver
        autoPtr<motionSolver> msPtr_;
        
        //- Is the valve really closed
        boolList isReallyClosed_;
       
        //- correct points motion after attach
        Switch correctPointsMotion_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        engineValveSliding(const engineValveSliding&);

        //- Disallow default bitwise assignment
        void operator=(const engineValveSliding&);

        //- Make layering modifiers live
        void makeLayersLive();

        //- Make Sliding interface modifiers live
        void makeSlidersLive();

        //- Is the sliding interface attached?
        bool attached() const;
        
        //- Prepare valve attach/detach
        void valveDetach();

        //- Prepare valve attach/detach
        void valveAttach();

        //- Prepare valve attach/detach
        void prepareValveDetach();

        //- Check if all patches exist, then calculate virtualPistonPosition, 
        //- pistonPosition and deckHeight for the first time
        void checkAndCalculate();
        
        //- Calculate the virtualPistonPosition, 
        void setVirtualPositions();
        
        //- Correct the vertical motion for the points below the "virtual piston", 
        void correctVerticalMotion();
        
        bool realDeformation() const;
                        
        bool inValve(const point& p, const label& i) const;
        bool inPiston(const point& p) const;

        bool isACylinderHeadFace(const labelList& cylHeadFaces, const label face);
                        
public:

    //- Runtime type information
    TypeName("engineValveSliding");


    // Constructors

        //- Construct from database
        explicit engineValveSliding(const IOobject& io);


    // Destructor - default


    // Member Functions

        inline const enginePiston& piston() const;
        inline const thoboisSlidingValveBank& valves() const;
        inline const scalar& deformSwitch() const;
        inline const scalar& valveTopTol() const;
        inline const scalar& pistonPosition() const;
        inline scalar& pistonPosition();
        inline const scalar& deckHeight() const;
        inline scalar& deckHeight();
        inline const scalarField& minValveZ() const;
        inline scalarField& maxValveZ();
        inline const scalar& poppetValveTol() const;
        inline boolList& isReallyClosed();
        inline const boolList& isReallyClosed() const;


        //- Return true for mesh deformation mode
        bool deformation() const
        {
            return
                engTime().thetaRevolution() > -deformSwitch_
             && engTime().thetaRevolution() < deformSwitch_;
        }
        
        //- Return number of valves
        label nValves() const
        {
            return valves_.size();
        }

        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers();

        //- Move and morph
        virtual bool update();
        
        //- Set boundary velocities
        void setBoundaryVelocity(volVectorField& U);
        
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "engineValveSlidingI.H"

#endif

// ************************************************************************* //
