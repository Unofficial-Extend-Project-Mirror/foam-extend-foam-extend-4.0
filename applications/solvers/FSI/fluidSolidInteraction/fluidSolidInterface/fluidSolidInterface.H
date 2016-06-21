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
    fluidSolidInterface

Description
    Fluid-structure interface class

Author
    Zeljko Tukovic, FSB Zagreb.  All rights reserved.

SourceFiles
    fluidSolidInterface.C

\*---------------------------------------------------------------------------*/

#ifndef fluidSolidInterface_H
#define fluidSolidInterface_H

#include "fluidSolver.H"
#include "solidSolver.H"
#include "IOdictionary.H"
#include "patchToPatchInterpolation.H"
#include "dynamicFvMesh.H"
#include "ggiInterpolation.H"
#include "movingWallVelocityFvPatchVectorField.H"

#include "RBFInterpolation.H"
#include "TPSFunction.H"
using namespace rbf;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class fluidSolidInterface Declaration
\*---------------------------------------------------------------------------*/

class fluidSolidInterface
:
    public IOdictionary
{
    // Private data

        //- Fluid (flow) mesh
        dynamicFvMesh& fluidMesh_;

        //- Fluid model
        autoPtr<fluidSolver> flow_;

        //- Solid (stress) mesh
        fvMesh& solidMesh_;

        //- Stress model
        autoPtr<solidSolver> stress_;

        //- Solid patch index
        label solidPatchIndex_;

        //- Solid face zone index
        label solidZoneIndex_;

        //- Fluid patch index
        label fluidPatchIndex_;

        //- Fluid face zone index
        label fluidZoneIndex_;

        //- Solid face zone current points
        mutable vectorField* currentSolidZonePointsPtr_;

        //- Solid face zone primitive patch in current configuration
        mutable PrimitivePatch<face, List, const pointField&>*
        currentSolidZonePatchPtr_;

        //- Fluid zone to solid zone interpolator
        mutable std::shared_ptr<RBFInterpolation> fluidToSolidPtr_;

//         //- Ggi zone-to-zone interpolation
//         mutable ggiZoneInterpolation* ggiFluidToSolidPtr_;

        //- Ggi zone-to-zone interpolation
        mutable ggiZoneInterpolation* ggiInterpolatorPtr_;

        //- Solid zone to fluid zone interpolator
        mutable std::shared_ptr<RBFInterpolation> solidToFluidPtr_;

        //- Coupling scheme
        word couplingScheme_;

        //- Relaxation factor
        scalar relaxationFactor_;

        //- Relaxation factor
        scalar aitkenRelaxationFactor_;

        //- Outer correction loop stoping tolerance
        scalar outerCorrTolerance_;

        //- Maximal number of outer correctors
        label nOuterCorr_;

        //- Fsi coupling on/off
        Switch coupled_;

        //- Predictor on/off
        Switch predictor_;

        //- Use rbf interpolation
        Switch rbfInterpolation_;

        //- Coupling reuse
        label couplingReuse_;

        //- Interface deformation limit
        scalar interfaceDeformationLimit_;

        //- Fluid zone point displacement
        vectorField fluidZonePointsDispl_;

        //- Fluid zone ref. point displacement
        vectorField fluidZonePointsDisplRef_;

        //- Fluid zone previous point displacement
        vectorField fluidZonePointsDisplPrev_;

        //- Solid zone point displacement
        vectorField solidZonePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField solidZonePointsDisplRef_;

        //- Solid zone point displacement
        vectorField interfacePointsDispl_;

        //- Solid zone ref. point displacement
        vectorField interfacePointsDisplPrev_;

        //- Solid zone pressure
        scalarField solidZonePressure_;

        //- Solid zone traction (pressure + vicous)
        vectorIOField solidZoneTraction_;

        //- Solid zone traction (pressure + vicous)
        vectorField solidZoneTractionPrev_;

        //- Solid zone traction (pressure + vicous)
        vectorIOField predictedSolidZoneTraction_;

        //- Current fsi residual
        vectorField residual_;

        //- Previous fsi residual
        vectorField residualPrev_;

        //- Maximal resudual norm
        scalar maxResidualNorm_;

        //- Maximal interface displacement norm
        scalar maxIntDisplNorm_;

        //- Outer corrector
        label outerCorr_;

//         //- Is it fluid domain pure Dirichlet (witout outlets)?
//         Switch closedFluidDomain_;

//         //- Reference pressure
//         scalar refPressure_;

//         //- Reference pressure
//         scalar refPressureIncrement_;

//         //- Current time index
//         label timeIndex_;

//         //- Copressibility
//         scalar compressibility_;

        //- Interpolator update frequency
        label interpolatorUpdateFrequency_;


        //- IQN-ILS coupling fields
        DynamicList<vectorField> fluidPatchPointsV_;
        DynamicList<vectorField> fluidPatchPointsW_;
        DynamicList<scalar> fluidPatchPointsT_;

        //- Accumulated fluid side interface displacement
        mutable vectorIOField* accumulatedFluidInterfaceDisplacementPtr_;

        //- Min edge length for interface points at fluid side
        mutable scalarField* minEdgeLengthPtr_;

    // Private Member Functions

        //- Calculate current solid zone points
        void calcCurrentSolidZonePoints() const;

        //- Calculate current solid zone primitive patch
        void calcCurrentSolidZonePatch() const;

        //- Calculate fluid to solid interpolator
        void calcFluidToSolidInterpolator() const;

//         //- Calculate fluid to solid ggi interpolator
//         void calcGgiFluidToSolidInterpolator() const;

        //- Calculate fluid to solid ggi interpolator
        void calcGgiInterpolator() const;

        //- Calculate fluid to solid interpolator
        void calcSolidToFluidInterpolator() const;

        //- Accumulated fluid interface displacement
        void calcAccumulatedFluidInterfaceDisplacement() const;

        //- Calculate minimal edge lengths
        void calcMinEdgeLength() const;

        //- Return accumulated interface displacement
        vectorIOField& accumulatedFluidInterfaceDisplacement();

        //- Return minimal edge length
        const scalarField& minEdgeLength() const;

        //- Disallow default bitwise copy construct
        fluidSolidInterface(const fluidSolidInterface&);

        //- Disallow default bitwise assignment
        void operator=(const fluidSolidInterface&);


protected:

    // Protected member functions

public:

    //- Runtime type information
    TypeName("fluidSolidInterface");


    // Declare run-time constructor selection table


    // Constructors

        //- Construct from components
        fluidSolidInterface
        (
            dynamicFvMesh& fluidMesh,
            fvMesh& solidMesh
        );


    // Selectors


    // Destructor

        virtual ~fluidSolidInterface();


    // Member Functions

        // Access

            //- Return fluid mesh
            const dynamicFvMesh& fluidMesh() const
            {
                return fluidMesh_;
            }

            //- Return solid mesh
            const fvMesh& solidMesh() const
            {
                return solidMesh_;
            }

            //- Return time
            const Time& runTime() const
            {
                return fluidMesh_.time();
            }

            //- Return fluid solver
            const fluidSolver& flow() const
            {
                return flow_();
            }

            //- Return fluid solver
            fluidSolver& flow()
            {
                return flow_();
            }

            //- Return solid solver
            const solidSolver& stress() const
            {
                return stress_();
            }

            //- Return solid solver
            solidSolver& stress()
            {
                return stress_();
            }

            //- Return solid patch index
            label solidPatchIndex() const
            {
                return solidPatchIndex_;
            }

            //- Return solid face zone index
            label solidZoneIndex() const
            {
                return solidZoneIndex_;
            }

            //- Return fluid patch index
            label fluidPatchIndex() const
            {
                return fluidPatchIndex_;
            }

            //- Return fluid face zone index
            label fluidZoneIndex() const
            {
                return fluidZoneIndex_;
            }

            //- Return current solid zone points
            const vectorField& currentSolidZonePoints() const;

            //- Return current solid zone patch
            const PrimitivePatch<face, List, const pointField&>&
            currentSolidZonePatch() const;

            //- Return fluid to solid interpolator
            const std::shared_ptr<RBFInterpolation>& fluidToSolid() const;

//             //- Return fluid to solid interpolator
//             const ggiZoneInterpolation& ggiFluidToSolid() const;

            //- Return fluid to solid interpolator
            const ggiZoneInterpolation& ggiInterpolator() const;

            //- Return fluid to solid interpolator
            const std::shared_ptr<RBFInterpolation>& solidToFluid() const;

            //- Return coupling scheme
            const word& couplingScheme() const
            {
                return couplingScheme_;
            }

            //- Return relaxation factor
            scalar relaxationFactor() const
            {
                return relaxationFactor_;
            }

            //- Return relaxation factor
            scalar& aitkenRelaxationFactor()
            {
                return aitkenRelaxationFactor_;
            }

            //- Return outer corrector loop tolerance
            scalar outerCorrTolerance() const
            {
                return outerCorrTolerance_;
            }

            //- Return max numter of outer correctors
            label nOuterCorr() const
            {
                return nOuterCorr_;
            }

            //- Is it fluid and structure coupled
            const Switch& coupled() const
            {
                return coupled_;
            }

            //- Is it fluid and structure coupled
            const Switch& predictor() const
            {
                return predictor_;
            }

            //- Is it fluid and structure coupled
            Switch& coupled()
            {
                return coupled_;
            }

            //- Is it fluid and structure coupled
            label couplingReuse() const
            {
                return couplingReuse_;
            }

            //- Return relaxation factor
            scalar interfaceDeformationLimit() const
            {
                return interfaceDeformationLimit_;
            }

            //- Return fluid zone point displacement
            vectorField& fluidZonePointsDispl()
            {
                return fluidZonePointsDispl_;
            }

            //- Return fluid zone ref. point displacement
            vectorField& fluidZonePointsDisplRef()
            {
                return fluidZonePointsDisplRef_;
            }

            //- Return fluid zone previous point displacement
            vectorField& fluidZonePointsDisplPrev()
            {
                return fluidZonePointsDisplPrev_;
            }

            //- Return solid zone point displacement
            vectorField& solidZonePointsDispl()
            {
                return solidZonePointsDispl_;
            }

            //- Return solid zone ref. point displacement
            vectorField& solidZonePointsDisplRef()
            {
                return solidZonePointsDisplRef_;
            }

            //- Return solid zone pressure
            scalarField& solidZonePressure()
            {
                return solidZonePressure_;
            }

            //- Return solid zone total traction
            vectorField& solidZoneTraction()
            {
                return solidZoneTraction_;
            }

            //- Return solid zone total traction
            vectorField& solidZoneTractionPrev()
            {
                return solidZoneTractionPrev_;
            }

            //- Return current fsi residual
            vectorField& residual()
            {
                return residual_;
            }

            //- Return previous fsi residual
            vectorField& residualPrev()
            {
                return residualPrev_;
            }

            //- Return current outer iteration
            label& outerCorr()
            {
                return outerCorr_;
            }

            //- Return current outer iteration
            const label& outerCorr() const
            {
                return outerCorr_;
            }

//             //- Is it fluid domain pure Dirichlet
//             Switch closedFluidDomain() const
//             {
//                 return closedFluidDomain_;
//             }

//             //- Ref. pressure used in case of closed fluid domain
//             scalar refPressure()
//             {
//                 return refPressure_;
//             }

//             //- Ref. pressure used in case of closed fluid domain
//             scalar refPressureIncrement()
//             {
//                 return refPressureIncrement_;
//             }

//             //- Compressibility used in case of closed fluid domain
//             scalar compressibility()
//             {
//                 return compressibility_;
//             }

        // Edit

            //- Initialize fields
            void initializeFields();

            //- Initialize fields
            void updateInterpolator();

            //- Calculate interface displacement
            void updateDisplacement();

            //- Calculate interface displacement
            scalar updateWeakDisplacement();

            //- Calculate interface displacement
            void updateDisplacementUsingIQNILS();

            //- Move fluid mesh
            void moveFluidMesh();

            //- Update interface force
            void updateForce();

            //- Update interface force
            void updateWeakForce();

            //- Update interface force
            void updateWeakTraction();

            //- Update interface force
            void predictAndUpdateForce();

            //- Solve solid
            void evolveStress();

            //- Update interface force
            scalar updateResidual();

//             //- Evolve the fluid solver
//             virtual void evolve() = 0;

//             //- Read dictionary
//             virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
