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

CLASS
    accordionEngineMesh

DESCRIPTION

*/
// ------------------------------------------------------------------------- //


#ifndef accordionEngineMesh_H
#define accordionEngineMesh_H

#include "engineTopoChangerMesh.H"
#include "enginePiston.H"
#include "motionSolver.H"
#include "polyPatchID.H"
#include "accordionValveBank.H"
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
        accordionValveBank valves_;

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
//        word pistonAuxPoints_;

        //- PistonFaceSet
        word headCellSetName_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        accordionEngineMesh(const accordionEngineMesh&);

        //- Disallow default bitwise assignment
        void operator=(const accordionEngineMesh&);

        //- Check if all patches exist, then calculate virtualPistonPosition,
        //- pistonPosition and deckHeight for the first time
        void checkAndCalculate();

        //- Calculate the virtualPistonPosition
        void setVirtualPositions();

        bool realDeformation() const;


public:

    //- Runtime type information
    TypeName("accordionEngineMesh");


    // Constructors

        //- Construct from database
        explicit accordionEngineMesh(const IOobject& io);


    // Destructor

        virtual ~accordionEngineMesh()
        {}


    // Member Functions

        inline const enginePiston& piston() const;
        inline const accordionValveBank& valves() const;
        inline const scalar& pistonPosition() const;
        inline scalar& pistonPosition();
        inline const scalar& virtualPistonPosition() const;
        inline scalar& virtualPistonPosition();
        inline const scalar& deckHeight() const;
        inline scalar& deckHeight();


        //- Return true for mesh deformation mode
        bool deformation() const
        {
            return true;
        }

        //- Return number of valves
        label nValves() const
        {
            return valves_.size();
        }

        //- Add valve and piston zones and modifiers
        void addMeshZones();

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
