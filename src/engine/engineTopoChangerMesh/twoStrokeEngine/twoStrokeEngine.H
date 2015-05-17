/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

\*---------------------------------------------------------------------------*/

#ifndef twoStrokeEngine_H
#define twoStrokeEngine_H

#include "engineTopoChangerMesh.H"
#include "enginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                      Class twoStrokeEngine Declaration
\*---------------------------------------------------------------------------*/

class twoStrokeEngine
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        enginePiston piston_;

        //- Scavenging ports patches
        wordList scavInPortPatches_;

        //- Scavenging cylinder patches
        wordList scavInCylPatches_;

        //- Name of the head pointSet
        word headPointsSetName_;

        //- Name of the head pointSet
        word headCellsSetName_;

        //- Name of the moving cell set
        word movingCellSetName_;


        //- Markup field for points.  Moving points marked with 1
        mutable scalarField* movingPointsMaskPtr_;

        //- Layering-to-deformation switch in crank-angle degrees
        scalar deformSwitch_;

        //- Piston Position
        scalar pistonPosition_;

        //- Virtual piston position (pistonPosition + offSet)
        scalar virtualPistonPosition_;

        //- deckHeight
        scalar deckHeight_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        twoStrokeEngine(const twoStrokeEngine&);

        //- Disallow default bitwise assignment
        void operator=(const twoStrokeEngine&);

        //- Make layering modifiers live
        void makeLayersLive();

        //- Make Sliding interface modifiers live
        void makeSlidersLive();

        //- Is the sliding interface attached?
        bool attached() const;

        //- Check if all patches exist, then calculate virtualPistonPosition,
        //- pistonPosition and deckHeight for the first time
        void checkAndCalculate();

        //- Calculate the virtualPistonPosition,
        void setVirtualPistonPosition();

        //- Correct the vertical motion for the points below the
        //  "virtual piston",
        void correctVerticalMotion();

        //- Apply the moving mask, the points of the ports do not move.
        void applyMovingMask();

        //- Calculate moving masks
        void calcMovingMasks() const;

        //- Return moving points mask
        const scalarField& movingPointsMask() const;

        bool realDeformation() const;

        void checkMotionFluxes();


public:

    //- Runtime type information
    TypeName("twoStrokeEngine");


    // Constructors

        //- Construct from database
        explicit twoStrokeEngine(const IOobject& io);


    // Destructor - default


    // Member Functions

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
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

        const enginePiston& piston() const
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

#endif

// ************************************************************************* //
