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

\*---------------------------------------------------------------------------*/

#ifndef simpleEngineTopoFvMesh_H
#define simpleEngineTopoFvMesh_H

#include "topoChangerFvMesh.H"
#include "engineTime.H"
#include "valveBank.H"
#include "simpleEnginePiston.H"
#include "motionSolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class simpleEngineTopoFvMesh Declaration
\*---------------------------------------------------------------------------*/

class simpleEngineTopoFvMesh
:
    public topoChangerFvMesh
{
    // Private data

        //- Engine database
        const engineTime& engineTime_;

        //- Engine valves
        valveBank valves_;

        //- Piston patch
        simpleEnginePiston piston_;

        //- Motion solver
        autoPtr<motionSolver> msPtr_;

        //- Layering-to-deformation switch in crank-angle degrees
        scalar deformSwitch_;

        //- Valve position tolerance.  Used in mask calculations
        scalar valvePosTol_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        simpleEngineTopoFvMesh(const simpleEngineTopoFvMesh&);

        //- Disallow default bitwise assignment
        void operator=(const simpleEngineTopoFvMesh&);


        //- Make layering modifiers live
        void makeLayersLive();

        //- Make sliding modifiers live
        void makeSlidersLive();

        //- Prepare valve attach/detach
        void prepareValveDetach();

        //- Return true if sliders are attached
        bool attached() const;

        //- Set valve and piston motion boundary conditions for deformation
        void setBoundaryMotion();

        //- Set valve and piston position for layering
        void setBoundaryPosition();


public:

    //- Runtime type information
    TypeName("simpleEngineTopoFvMesh");


    // Constructors

        //- Construct from database
        explicit simpleEngineTopoFvMesh(const IOobject& io);


    // Destructor

        virtual ~simpleEngineTopoFvMesh()
        {}


    // Member Functions

        //- Return valve bank
        const valveBank& valves() const
        {
            return valves_;
        }

        //- Return number of valves
        label nValves() const
        {
            return valves().size();
        }

        const simpleEnginePiston& piston() const
        {
            return piston_;
        }

        //- Return true for mesh deformation mode
        bool deformation() const
        {
            return
                engineTime_.thetaRevolution() > -deformSwitch_
             && engineTime_.thetaRevolution() < deformSwitch_;
        }

        //- Add valve and piston zones and modifiers
        void addZonesAndModifiers();

        //- Update the mesh for both mesh motion and topology change
        virtual bool update();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
