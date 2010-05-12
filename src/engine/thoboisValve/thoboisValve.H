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

Class
    thoboisValve

Description

SourceFiles
    thoboisValve.C

\*---------------------------------------------------------------------------*/

#ifndef thoboisValve_H
#define thoboisValve_H

#include "word.H"
#include "coordinateSystem.H"
#include "polyPatchID.H"
#include "graph.H"

#include "faceSet.H"
#include "pointSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Class forward declarations
class polyMesh;
class engineTime;

/*---------------------------------------------------------------------------*\
                          Class thoboisValve Declaration
\*---------------------------------------------------------------------------*/

class thoboisValve
{
    // Private data

        //- Name of valve
        word name_;

        //- Reference to engine mesh
        const polyMesh& mesh_;

        //- Reference to engine database
        const engineTime& engineDB_;

        //- Coordinate system
        autoPtr<coordinateSystem> csPtr_;


        // Patch and zone names

            //- Valve bottom patch
            polyPatchID bottomPatch_;

            //- Valve poppet patch
            polyPatchID poppetPatch_;

            //- Valve poppet patch
            polyPatchID sidePatch_;

            //- Valve stem patch
            polyPatchID stemPatch_;

            //- Valve detach in cylinder patch
            polyPatchID detachInCylinderPatch_;

            //- Valve detach in port patch
            polyPatchID detachInPortPatch_;

            //- Faces to detach
            //labelList detachFaces_;
	    word detachFacesName_;
           
        // Valve lift data

            //- Valve lift profile
            graph liftProfile_;

            //- Lift curve start angle
            scalar liftProfileStart_;

            //- Lift curve end angle
            scalar liftProfileEnd_;

            //- Minimum valve lift.  On this lift the valve is considered closed
            const scalar minLift_;

            //- Valve diameter
            const scalar diameter_;

        // Mesh motion data
            
            //- Name of the pointSet for the static points
            word staticPointsName_;
            
            //- Name of the pointSet for the moving points
            word movingPointsName_;

            //- Name of the pointSet for the moving internal points
            word movingInternalPointsName_;
            
            //- Name of the cellSet for the static cells
            word staticCellsName_;
            
            //- Name of the cellSet for the moving cells
            word movingCellsName_;
        
        
    // Private Member Functions

        //- Disallow default bitwise copy construct
        thoboisValve(const thoboisValve&);

        //- Disallow default bitwise assignment
        void operator=(const thoboisValve&);


        //- Adjust crank angle to drop within the limits of the lift profile
        scalar adjustCrankAngle(const scalar theta) const;

public:

    // Static data members


    // Constructors

        //- Construct from components
        thoboisValve
        (
            const word& name,
            const polyMesh& mesh,
            const autoPtr<coordinateSystem>& valveCS,
            const word& bottomPatchName,
            const word& poppetPatchName,
            const word& sidePatchName,
            const word& stemPatchName,
            const word& detachInCylinderPatchName,
            const word& detachInPortPatchName,
            const word& detachFacesName,
            const graph& liftProfile,
            const scalar minLift,
            const scalar diameter,
            const word& staticPointsName,
            const word& movingPointsName,
            const word& movingInternalPointsName,
            const word& staticCellsName,
            const word& movingCellsName
        );

        //- Construct from dictionary
        thoboisValve
        (
            const word& name,
            const polyMesh& mesh,
            const dictionary& dict
        );


    // Destructor - default


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return coordinate system
        const coordinateSystem& cs() const
        {
            return csPtr_();
        }

        //- Return lift profile
        const graph& liftProfile() const
        {
            return liftProfile_;
        }

        //- Return valve diameter
        scalar diameter() const
        {
            return diameter_;
        }


        // Valve patches

            //- Return ID of bottom patch
            const polyPatchID& bottomPatchID() const
            {
                return bottomPatch_;
            }

            //- Return ID of poppet patch
            const polyPatchID& poppetPatchID() const
            {
                return poppetPatch_;
            }

            //- Return ID of side patch
            const polyPatchID& sidePatchID() const
            {
                return sidePatch_;
            }

            //- Return ID of stem patch
            const polyPatchID& stemPatchID() const
            {
                return stemPatch_;
            }

            //- Return ID of detach in cylinder patch
            const polyPatchID& detachInCylinderPatchID() const
            {
                return detachInCylinderPatch_;
            }

            //- Return ID of detach in port patch
            const polyPatchID& detachInPortPatchID() const
            {
                return detachInPortPatch_;
            }

            //- Return face labels of detach curtain
            const word& detachFacesName() const
            {
                return detachFacesName_;
            }
            
            const word& staticPointsName() const
            {
                return staticPointsName_;
            }

            const word& movingPointsName() const
            {
                return movingPointsName_;
            }

            const word& movingInternalPointsName() const
            {
                return movingInternalPointsName_;
            }

            const word& staticCellsName() const
            {
                return staticCellsName_;
            }

            const word& movingCellsName() const
            {
                return movingCellsName_;
            }

        // Valve position and velocity

            //- Return valve lift given crank angle in degrees
            scalar lift(const scalar theta) const;

            //- Is the valve open?
            bool isOpen() const;

            //- Return current lift
            scalar curLift() const;

            //- Return valve velocity for current time-step
            scalar curVelocity() const;

            //- Return list of active patch labels for the valve head
            //  (stem is excluded)
            labelList movingPatchIDs() const;


        //- Write dictionary
        void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
