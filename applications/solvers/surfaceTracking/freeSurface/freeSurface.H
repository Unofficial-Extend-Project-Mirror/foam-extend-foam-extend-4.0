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
    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

Class
    freeSurface

Description
    Implementation of Free surface (interface) tracking method 
    based on moving mesh.

SourceFiles
    freeSurface.C
    freeSurfacePointDisplacement.C
    makeFreeSurfaceData.C

\*---------------------------------------------------------------------------*/

#ifndef FreeSurface_H
#define FreeSurface_H


#include "fvCFD.H"
#include "IOdictionary.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "faCFD.H"
#include "patchToPatchInterpolation.H"
#include "tetDecompositionMotionSolver.H"
#include "tetPolyPatchInterpolation.H"

#include "surfactantProperties.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class freeSurface Declaration
\*---------------------------------------------------------------------------*/

class freeSurface
:
    public IOdictionary
{
    // Private data
    
        //- Reference to fvMesh 
        fvMesh& mesh_;

        //- Reference to density field
        const volScalarField& rho_;

        //- Reference to velocity field
        volVectorField& U_;

        //- Reference to pressure field
        volScalarField& p_;

        //- Reference to fluid flux field
        const surfaceScalarField& phi_;

        //- Current time index for free-surface
        label curTimeIndex_;

        //- Are there exist two fluids in the free-surface model
        Switch twoFluids_;

        //- Is it free-surface points displacement direction
        //  parallel with free-surface point normals
        Switch normalMotionDir_;

        //- Is it necessary to smoot the free surface
        Switch freeSurfaceSmoothing_;

        //- Is it free-surface clean of surfactants
        Switch cleanInterface_;

        //- ID of the fvMesh patch which represents 
        //  A side of the free-surface
        label aPatchID_; 

        //- ID of the fvMesh patch which represents 
        //  B side of the free-surface
        label bPatchID_; 

        //- Viscosity of fluid A
        dimensionedScalar muFluidA_;
    
        //- Viscosity of fluid B
        dimensionedScalar muFluidB_;

        //- Density of fluid A
        dimensionedScalar rhoFluidA_;
    
        //- Density of fluid B
        dimensionedScalar rhoFluidB_;

        //- Reference to gravitational acceleration
        dimensionedVector g_;

        //- Surface tension for the clean free-surface
        dimensionedScalar cleanInterfaceSurfTension_;

        //- Free surface patches which do not move
        wordList fixedFreeSurfacePatches_;

    // Demand-driven data

        //- Patch A motion interpolation
        mutable tetPolyPatchInterpolation* interpolatorAPtr_;

        //- Patch B motion interpolation
        mutable tetPolyPatchInterpolation* interpolatorBPtr_;

        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface 
        //  (interface) patches
        mutable patchToPatchInterpolation* interpolatorABPtr_;
    
        //- Patch to patch interpolation object which deals with
        //  interpolation of flow variables between two free-surface 
        //  (interface) patches
        mutable patchToPatchInterpolation* interpolatorBAPtr_;

        //- Points which are attached to the free-surface A side faces
        //  and which defines the free-surface shape
        mutable vectorIOField* controlPointsPtr_;

        //- Field which additionally determines 
        //  the motion of free-surface points
        mutable scalarField* motionPointsMaskPtr_;

        //- Displacement direction of free-surface points
        mutable vectorField* pointsDisplacementDirPtr_;

        //- Displacement direction of free-surface control points
        mutable vectorField* facesDisplacementDirPtr_;

        //- Total displacement of free-surface points 
        //  in one time step
        mutable vectorIOField* totalDisplacementPtr_;

        //- Finite area mesh which consists
        //  free-surface side A faces
        mutable faMesh* aMeshPtr_;

        //- Reference to mesh motion solver
        mutable tetDecompositionMotionSolver* mSolverPtr_;

        //- Free-surface velocity field
        mutable areaVectorField* UsPtr_;

        //- Free-surface fluid flux
        mutable edgeScalarField* phisPtr_;

        //- Free-surface surfactant concetration
        mutable areaScalarField* surfactConcPtr_;

        //- Surface tension field
        mutable areaScalarField* surfaceTensionPtr_;

        //- Surfactant properties
        mutable surfactantProperties* surfactantPtr_;

        //- Fluid indicator
        mutable volScalarField* fluidIndicatorPtr_;

    // Private Member Functions

        // Make demand-driven data

            void makeInterpolators();
            void makeControlPoints();
            void makeMotionPointsMask();
            void makeDirections();
            void makeTotalDisplacement();
            void readTotalDisplacement();
            void makeFaMesh();
            void makeMeshMotionSolver();
            void makeUs();
            void makePhis();
            void makeSurfactConc();
            void makeSurfaceTension();
            void makeSurfactant();
            void makeFluidIndicator();
            void makePatchPointInterpolators();


        //- Clear all demand-driven data
        void clearOut();

        //- Disallow default bitwise copy construct
        freeSurface(const freeSurface&);

        //- Disallow default bitwise assignment
        void operator=(const freeSurface&);

        //- Move control ponts for deltaH and calculate
        //  free surface points displacement
        //  for given new control points position
        tmp<vectorField> pointDisplacement(const scalarField& deltaH);

public:

    // Declare name of the class and it's debug switch
    ClassName("freeSurface");


    // Constructors

    freeSurface
    (
        fvMesh& m,
        const volScalarField& rho,
        volVectorField& Ub, 
        volScalarField& Pb, 
        const surfaceScalarField& phi
    );


    // Destructor

    ~freeSurface();  


    // Member Functions


        //- Return reference to fvMesh
        const fvMesh& mesh() const
        {
            return mesh_;
        }

        //- Return reference to fvMesh
        fvMesh& mesh()
        {
            return mesh_;
        }

        //- Return reference to time
        const Time& time() const
        {
            return U_.mesh().time();
        }

        //- Return reference to velocity field
        const volScalarField& rho() const
        {
            return rho_;
        }

        //- Return reference to velocity field
        volVectorField& U()
        {
            return U_;
        }

        //- Return reference to pressure field
        volScalarField& p()
        {
            return p_;
        }

        //- Return reference to flux field
        const surfaceScalarField& phi() const
        {
            return phi_;
        }

        //- Are there exist two fluid
        //  in the free-surface model
        bool twoFluids() const
        {
            return twoFluids_;
        }

        //- Motion direction swithc
        bool normalMotionDir() const
        {
            return normalMotionDir_;
        }

        //- Is it necessary to consider the normal shear stress 
        //  at the free-surface
        bool freeSurfaceSmoothing() const
        {
            return freeSurfaceSmoothing_;
        }

        //- Is it interface clean
        bool cleanInterface() const
        {
            return cleanInterface_;
        }

        //- ID of patch which represent A side
        //  of free-surface
        const label& aPatchID() const
        {
            return aPatchID_;
        }

        //- ID of patch which represent B side
        //  of free-surface
        const label& bPatchID() const
        {
            return bPatchID_;
        }

        //- Viscosity of fluid A
        const dimensionedScalar& muFluidA() const
        {
            return muFluidA_;
        }

        //- Viscosity of fluid B
        const dimensionedScalar& muFluidB() const
        {
            return muFluidB_;
        }

        //- Density of fluid A
        const dimensionedScalar& rhoFluidA() const
        {
            return rhoFluidA_;
        }

        //- Densiti of fluid B
        const dimensionedScalar& rhoFluidB() const
        {
            return rhoFluidB_;
        }

        //- Gravity
        const dimensionedVector& g() const
        {
            return g_;
        }

        //- Surface tension for clean interface
        const dimensionedScalar& cleanInterfaceSurfTension() const
        {
            return cleanInterfaceSurfTension_;
        }

        //- Return reference to interpolator for point motion, side A
        const tetPolyPatchInterpolation& interpolatorA();

        //- Return reference to interpolator for point motion, side B
        const tetPolyPatchInterpolation& interpolatorB();

        //- Return reference to interpolator from free surface side A to
        //  free surface side B
        const patchToPatchInterpolation& interpolatorAB();

        //- Return reference to interpolator from free surface side B to
        //  free surface side A
        const patchToPatchInterpolation& interpolatorBA();

        //- Return control points
        vectorField& controlPoints();

        //- Return reference to motion points mask field
        scalarField& motionPointsMask();

        //- Return reference to point displacement direction field
        vectorField& pointsDisplacementDir();

        //- Return reference to control points displacement direction field
        vectorField& facesDisplacementDir();

        //- Return total points displacement
        vectorField& totalDisplacement();

        //- Return reference to finite area mesh
        faMesh& aMesh();

        //- Return reference to motion solver
        tetDecompositionMotionSolver& meshMotionSolver();
    
        //- Return free-surface velocity field
        areaVectorField& Us();

        //- Return free-surface fluid flux field
        edgeScalarField& Phis();

        //- Return free-surface surfactant concentration field
        areaScalarField& surfactantConcentration();

        //- Return surface tension field
        areaScalarField& surfaceTension();

        //- Return surfactant properties
        const surfactantProperties& surfactant();

        //- Return referenct to fluid indicator
        const volScalarField& fluidIndicator();

        //- Update control points end displacement directions
        void updateDisplacementDirections();

        //- Move only free-surface points
        bool movePoints();

        //- Move only free-surface points for uniform displacement 
        bool movePoints(scalar displ);

        //- Move mesh points for old time free-surface points displacement
        bool moveMeshPointsForOldFreeSurfDisplacement();

        //- Move mesh points
        bool moveMeshPoints();

        //- Move mesh points for given displacement filed
        bool moveMeshPoints(const scalarField& delta);

        //- Correct free-surface boundary conditions
        void correctBoundaryConditions();
    
        //- Correct free-surface pressure
        void correctPressure();

        //- Correct free-surface velocity
        void correctVelocity();

        //- Correct free-surface velocity gradient
        void correctVelocityGradient();

        //- Correct free-surface velocity
        void correctSurfaceFlux();

        //- Correct free-surface surfactant concentration 
        void correctSurfactantConcentration();

        //- Calculate total pressure force
        vector totalPressureForce();

        //- Calculate total viscous force
        vector totalViscousForce();

        //- Calculate free-surface undulation indicator
        tmp<scalarField> undulationIndicator();

        //- Smooth free-surface
        void smooth();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
