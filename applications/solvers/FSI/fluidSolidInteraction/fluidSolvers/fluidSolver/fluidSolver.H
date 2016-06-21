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
    fluidSolver

Description
    Virtual base class for fluid solvers

Author
    Zeljko Tukovic, Fsb Zagreb.  All rights reserved.

SourceFiles
    fluidSolver.C
    newFluidSolver.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSolver_H
#define fluidSolver_H

#include "fvMesh.H"
#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidSolver Declaration
\*---------------------------------------------------------------------------*/

class fluidSolver
:
    public IOdictionary
{
    // Private data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Fluid properties dictionary
        dictionary fluidProperties_;

        //- Global face zones
        mutable labelList* globalFaceZonesPtr_;

        //- Global to local face zone poin maping
        mutable labelListList* globalToLocalFaceZonePointMapPtr_;

    // Private Member Functions

        //- Find global face zones
        void calcGlobalFaceZones() const;

        //- Calculate global to local point addressing
        void calcGlobalToLocalFaceZonePointMap() const;

        //- Disallow default bitwise copy construct
        fluidSolver(const fluidSolver&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSolver&);


protected:

    // Protected member functions

        //- Return flow properties dictionary
        const dictionary& fluidProperties() const
        {
            return fluidProperties_;
        }


public:

    //- Runtime type information
    TypeName("fluidSolver");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            fluidSolver,
            dictionary,
            (
                const fvMesh& mesh
            ),
            (mesh)
        );


    // Constructors

        //- Construct from components
        fluidSolver
        (
            const word& type,
            const fvMesh& mesh
        );


    // Selectors

        //- Select constructed from mesh
        static autoPtr<fluidSolver> New(const fvMesh& mesh);


    // Destructor

        virtual ~fluidSolver();


    // Member Functions

        // Access

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

            //- Return velocity field
            virtual const volVectorField& U() const = 0;

            //- Return pressure field
            virtual const volScalarField& p() const = 0;

            //- Patch viscous force (N/m2)
            virtual tmp<vectorField> patchViscousForce
            (
                const label patchID
            ) const = 0;

            //- Patch pressure force (N/m2)
            virtual tmp<scalarField> patchPressureForce
            (
                const label patchID
            ) const = 0;

            //- Face zone viscous force (N/m2)
            virtual tmp<vectorField> faceZoneViscousForce
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone pressure force (N/m2)
            virtual tmp<scalarField> faceZonePressureForce
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Face zone velocity (m/s)
            virtual tmp<vectorField> faceZoneVelocity
            (
                const label zoneID,
                const label patchID
            ) const;

            //- Face zone effective dynamic viscosity
            virtual tmp<scalarField> faceZoneMuEff
            (
                const label zoneID,
                const label patchID
            ) const = 0;

            //- Get global face zones list
            const labelList& globalFaceZones() const;

            //- Get global to local face zone point addressing
            const labelListList& globalToLocalFaceZonePointMap() const;

        // Edit

            //- Evolve the fluid solver
            virtual void evolve() = 0;

            //- Update fields
            virtual void updateFields()
            {}

            //- Correct velocity at specified patch
            virtual void correctVelocity(label patchIndex)
            {}

            //- Read fluidProperties dictionary
            virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
