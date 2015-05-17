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

Class

Description
    Engine mesh class with dynamic layering on the piston patch.

Author
    Tommaso Lucchini, Politecnico di Milano.

SourceFiles
    addpistonLayerModifiers.C
    pistonLayer.C

\*---------------------------------------------------------------------------*/

#ifndef pistonLayer_H
#define pistonLayer_H

#include "engineTopoChangerMesh.H"
#include "enginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "polyMesh.H"
#include "volPointInterpolation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class pistonLayer Declaration
\*---------------------------------------------------------------------------*/

class pistonLayer
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        enginePiston piston_;

        //- Layering-to-deformation switch in crank-angle degrees
        scalar deformSwitch_;

        //- name of the head pointSet
        word headPointsSetName_;

        //- Piston Position
        scalar pistonPosition_;

        //- Virtual piston position (pistonPosition + offSet)
        scalar virtualPistonPosition_;

        //- deckHeight
        scalar deckHeight_;

        //- deckHeight
        Switch scalePoints_;

        //- deckHeight
        Switch movePointsBelowPiston_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        pistonLayer(const pistonLayer&);

        //- Disallow default bitwise assignment
        void operator=(const pistonLayer&);

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
    TypeName("pistonLayer");


    // Constructors

        //- Construct from database
        explicit pistonLayer(const IOobject& io);


    // Destructor

        virtual ~pistonLayer();

    // Member Functions


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
