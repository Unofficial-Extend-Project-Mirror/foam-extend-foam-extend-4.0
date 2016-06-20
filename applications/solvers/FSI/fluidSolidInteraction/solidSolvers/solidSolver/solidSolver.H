/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    solidSolver

Description
    Virtual base class for stress analysis models

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    solidSolver.C
    newSolidSolver.C

\*---------------------------------------------------------------------------*/

#ifndef solidSolver_H
#define solidSolver_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "constitutiveModel.H"
#include "materialInterface.H"



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class solidSolver Declaration
\*---------------------------------------------------------------------------*/

class solidSolver
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Flow properties dictionary
        dictionary solidProperties_;

        //- Body force
        volVectorField bodyForce_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone poin maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

        //- Number of outter correctors in last
        //  solution of solid solver
        label nCorrectors_;

        //- Use interpolated face gradient of displacement
        Switch interpFaceGrad_;

        //- Use least squares method to calculate cell
        //  gradient of displacement
        Switch lsCellGrad_;

    // Private Member Functions

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Find global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

        //- Disallow default bitwise copy construct
        solidSolver(const solidSolver&);

        //- Disallow default bitwise assignment
        void operator=(const solidSolver&);


protected:

    // Protected member functions

        //- Return flow properties dictionary
        dictionary& solidProperties()
        {
            return solidProperties_;
        }

public:

    //- Runtime type information
    TypeName("solidSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            solidSolver,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        solidSolver
        (
            const word& type,
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<solidSolver> New(const fvMesh& mesh);


    // Destructor

        virtual ~solidSolver();


    // Member Functions

        // Access

            //- Return face grad switch
            Switch interpFaceGrad() const
            {
                return interpFaceGrad_;
            }

            //- Return cell grad switch
            Switch lsCellGrad() const
            {
                return lsCellGrad_;
            }

            //- Return mesh
            const fvMesh& mesh() const
            {
                return mesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return mesh_.time();
            }

            //- Return flow properties dictionary
            const dictionary& solidProperties() const
            {
                return solidProperties_;
            }

            //- Thermal stress
            bool thermalStress() const;

            //- Return body force
            const volVectorField& bodyForce() const
            {
                return bodyForce_;
            }

            //- Return body force
            volVectorField& bodyForce()
            {
                return bodyForce_;
            }

            //- Constitutive model
            virtual const constitutiveModel& rheology() const = 0;

            //- Material interface
            virtual const autoPtr<materialInterface>& interface() const = 0;

            //- Material interface
            virtual autoPtr<materialInterface>& interface()  = 0;

            //- Return point displacement field
            virtual const pointVectorField& pointD() const = 0;

            //- Return point velocity field
            virtual vector pointU(label pointID) const = 0;

            //- Patch point displacement increment
            virtual tmp<vectorField> patchPointDisplacementIncrement
            (
                const label patchID
            ) const = 0;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacementIncrement
            (
                const label zoneID
            ) const = 0;

            //- Patch point displacement increment
            virtual tmp<vectorField> patchPointDisplacement
            (
                const label patchID
            ) const = 0;

            //- Face zone point displacement increment
            virtual tmp<vectorField> faceZonePointDisplacement
            (
                const label zoneID
            ) const = 0;

            //- Patch face acceleration
            virtual tmp<vectorField> patchFaceAcceleration
            (
                const label patchID
            ) const;

            //- Face zone velocity
            virtual tmp<vectorField> faceZoneAcceleration
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone velocity
            virtual tmp<vectorField> faceZoneVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone velocity gradient
            virtual tmp<tensorField> faceZoneSurfaceGradientOfVelocity
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone current normal
            virtual tmp<vectorField> faceZoneNormal
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone primitive patch in current configuration
            virtual tmp<vectorField>
            currentFaceZonePoints(const label zoneID) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

        // Edit

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const vectorField& traction
            ) = 0;

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID,
                const scalarField& pressure
            ) = 0;

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& faceZoneTraction
            ) = 0;

            //- Set traction at specified patch
            virtual void setPressure
            (
                const label patchID,
                const label zoneID,
                const scalarField& faceZonePressure
            ) = 0;

            //- Set traction at specified patch
            virtual void setVelocityAndTraction
            (
                const label patchID,
                const vectorField& traction,
                const vectorField& velocity,
                const vectorField& normal
            ){}

            //- Set traction at specified patch
            virtual void setVelocityAndTraction
            (
                const label patchID,
                const label zoneID,
                const vectorField& faceZoneTraction,
                const vectorField& faceZoneVelocity,
                const vectorField& faceZoneNormal
            ){}

            //- Predict traction at specified zone
            virtual tmp<vectorField> predictTraction
            (
                const label patchID,
                const label zoneID
            ) = 0;

            //- Predict pressure at specified zone
            virtual tmp<scalarField> predictPressure
            (
                const label patchID,
                const label zoneID
            ) = 0;

            //- Evolve the solid solver
            virtual bool evolve() = 0;

            //- Predict the solid solver
            virtual void predict() = 0;

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                label nCorrectors = 0
            ) = 0;

            //- Smooth solution
            virtual scalar smooth
            (
                const volVectorField& source,
                const volVectorField& refSolution,
                label nCorrectors = 0
            ) = 0;

            //- Return reference to current solution
            virtual const volVectorField& solution() const = 0;

            //- Return reference to current solution field
            virtual volVectorField& solution() = 0;

            //- Initialise solution
            virtual void initialise
            (
                const volVectorField& sol,
                bool consistenBoundaryField = false
            ) = 0;

            //- Return current residual
            virtual tmp<volVectorField> residual
            (
                const volVectorField& source
            ) = 0;


            //- Return number of outter correctors in last evaluation
            //  of solid solver
            label& nCorrectors()
            {
                return nCorrectors_;
            }

            //- Return number of outter correctors in last evaluation
            //  of solid solver
            const label& nCorrectors() const
            {
                return nCorrectors_;
            }

            //- Update total fields
            virtual void updateTotalFields()
            {}

            //- Update fields
            virtual void updateFields()
            {}

            //- Get current traction force
            virtual tmp<surfaceVectorField> traction() const = 0;

            //- Get current hydrosatic pressure
            virtual tmp<volScalarField> hydPressure() const;

            //- Read fluidProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
