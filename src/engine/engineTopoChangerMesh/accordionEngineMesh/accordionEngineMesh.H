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

CLASS
    accordionEngineMesh

DESCRIPTION
    Engine mesh class for canted valves engine. Layers are added/removed on the
    piston, while mesh is deformed close to the valves. Attach/detach boundary
    is used. There are no sliding interfaces.

*/
// ------------------------------------------------------------------------- //

#ifndef accordionEngineMesh_H
#define accordionEngineMesh_H

#include "engineTopoChangerMesh.H"
#include "enginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "thoboisValveBank.H"
#include "twoDPointCorrector.H"

#include "faceSet.H"
#include "pointSet.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations

/*---------------------------------------------------------------------------*\
                           Class lateralValves Declaration
\*---------------------------------------------------------------------------*/
class twoDPointCorrector;

class accordionEngineMesh
:
    public engineTopoChangerMesh
{
    // Private data

        //- Piston patch
        enginePiston piston_;
        
        //- Lateral valves
        thoboisValveBank valves_;
                         
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
                               
        //- Motion solver
        autoPtr<motionSolver> msPtr_;

        //- Cylinder head name
        word cylinderHeadName_;               
        
        //- Cylinder head name
        word linerName_;               
        
        //- PistonFaceSet
        word pistonAuxPoints_;
        
	Switch moveDetach_;
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        accordionEngineMesh(const accordionEngineMesh&);

        //- Disallow default bitwise assignment
        void operator=(const accordionEngineMesh&);

        //- Make layering modifiers live
        void makeLayersLive();

        //- Make attach-detach live
        void makeDetachLive();
       
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
                
        bool realDeformation() const;
                                                
public:

    //- Runtime type information
    TypeName("accordionEngineMesh");


    // Constructors

        //- Construct from database
        explicit accordionEngineMesh(const IOobject& io);


    // Destructor - default


    // Member Functions

        inline const enginePiston& piston() const;
        inline const thoboisValveBank& valves() const;
        inline const scalar& deformSwitch() const;
        inline const scalar& delta() const;
        inline const scalar& offSet() const;
        inline const scalar& pistonPosition() const;
        inline scalar& pistonPosition();
        inline const scalar& virtualPistonPosition() const;
        inline scalar& virtualPistonPosition();
        inline const scalar& deckHeight() const;
        inline scalar& deckHeight();


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

        //- Generate the mesh from a series of square blocks
        void generateMesh(){};
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "accordionEngineMeshI.H"

#endif

// ************************************************************************* //
