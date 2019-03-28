/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::solutionControl

Description
    Base class for solution control classes.
    The class also provides additional member functions for calculation
    of time and under-relaxation consistent flux and velocity.

\*---------------------------------------------------------------------------*/

#ifndef solutionControl_H
#define solutionControl_H

#include "fvMesh.H"
#include "surfaceMesh.H"
#include "fvsPatchField.H"
#include "fvMatrices.H"
#include "wordRe.H"
#include "MRFZones.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class solutionControl Declaration
\*---------------------------------------------------------------------------*/

class solutionControl
:
    public IOobject
{
public:

    struct fieldData
    {
        wordRe name;
        scalar absTol;
        scalar relTol;
        scalar initialResidual;
    };


protected:

    // Protected data

        //- Reference to the mesh database
        fvMesh& mesh_;

        //- List of residual data per field
        List<fieldData> residualControl_;

        //- The dictionary name, e.g. SIMPLE, PIMPLE
        const word algorithmName_;


        // Solution controls

            //- Maximum number of non-orthogonal correctors
            label nNonOrthCorr_;

            //- Flag to indicate to solve for momentum
            bool momentumPredictor_;

            //- Flag to indicate to solve using transonic algorithm
            bool transonic_;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            bool consistent_;


        // Evolution

            //- Current corrector loop index
            label corr_;

            //- Current non-orthogonal corrector loop index
            label corrNonOrtho_;


        // Fields necessary for time and under-relaxation consistency

            //- A coeff (A^~) arising from consistency formulation. Note: we can
            //  have multiple pressure/velocity systems, hence the PtrList
            mutable PtrList<surfaceScalarField> aCoeffPtrs_;

            //- Face velocity needed for consistent formulation. Note: we can
            //  have multiple pressure/velocity systems, hence the PtrList
            mutable PtrList<surfaceVectorField> faceUPtrs_;

            //- Hash Table containing indices of PtrLists for given names
            mutable HashTable<label> indices_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        virtual void read(const bool absTolOnly);

        //- Read controls from fvSolution dictionary
        virtual void read();

        //- Return index of field in residualControl_ if present
        virtual label applyToField
        (
            const word& fieldName,
            const bool useRegEx = true
        ) const;

        //- Return true if all convergence checks are satisfied
        virtual bool criteriaSatisfied() = 0;

        //- Store previous iteration fields
        virtual void storePrevIterFields() const;

        //- Store previous iteration field for vol<Type>Fields
        template<class Type>
        void storePrevIter() const;

        template<class Type>
        void maxTypeResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& firstRes,
            scalar& lastRes
        ) const;

        scalar maxResidual
        (
            const word& fieldName,
            ITstream& data,
            scalar& lastRes
        ) const;


        // Time and under-relaxation consistency helper functions

            //- Get relaxation factor for velocity field. Overriden in
            //  pimpleControl since different relaxation factor may be used for
            //  final iteration.
            virtual const dimensionedScalar relaxFactor
            (
                const volVectorField& U
            ) const;

            //- Add consistent flux contribution arising from time derivative
            //  term. Note: aCoeff is parameter for clarity
            void addDdtFluxContribution
            (
                surfaceScalarField& phi,
                surfaceScalarField& aCoeff,
                const surfaceVectorField& faceU,
                const volVectorField& U,
                const surfaceScalarField& rAUf,
                const fvVectorMatrix& ddtUEqn
            ) const;

            //- Add consistent flux contribution arising from the
            //  under-relaxation. Note: aCoeff is parameter for clarity
            void addUnderRelaxationFluxContribution
            (
                surfaceScalarField& phi,
                surfaceScalarField& aCoeff,
                const volVectorField& U
            ) const;

            //- Correct flux at the boundaries
            void correctBoundaryFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U
            ) const;

            void correctBoundaryFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U,
                const surfaceScalarField& meshPhi
            ) const;


private:

    // Private member functions

        //- Disallow default bitwise copy construct
        solutionControl(const solutionControl&);

        //- Disallow default bitwise assignment
        void operator=(const solutionControl&);


public:


    // Static Data Members

        //- Run-time type information
        TypeName("solutionControl");


    // Constructors

        //- Construct from mesh
        solutionControl(fvMesh& mesh, const word& algorithmName);


    //- Destructor
    virtual ~solutionControl();


    // Member Functions

        // Access

            //- Return the solution dictionary
            inline const dictionary& dict() const;

            //- Current corrector loop index
            inline label corr() const;

            //- Current non-orthogonal corrector index
            inline label corrNonOrtho() const;


        // Solution control

            //- Maximum number of non-orthogonal correctors
            inline label nNonOrthCorr() const;

            //- Helper function to identify final non-orthogonal iteration
            inline bool finalNonOrthogonalIter() const;

            //- Flag to indicate to solve for momentum
            inline bool momentumPredictor() const;

            //- Flag to indicate to solve using transonic algorithm
            inline bool transonic() const;

            //- Flag to indicate to relax pressure using the
            //  "consistent" approach of SIMPLEC
            inline bool consistent() const;


        // Time and under-relaxation consistency.
        // Note: argument matching U parameter needs to be equal U = H/A, where
        // H and A come from convection-difussion equation only (without time
        // derivative term)

            //- Calculate consistent flux for transient solvers (before pressure
            //  equation).
            void calcTransientConsistentFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U,
                const volScalarField& rAU,
                const fvVectorMatrix& ddtUEqn
            ) const;

            //- Calculate consistent flux for steady state solvers (before
            //  pressure equation).
            void calcSteadyConsistentFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U
            ) const;

            void calcSteadyConsistentFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U,
                const surfaceScalarField& meshPhi
            ) const;

            void calcSteadyMRFConsistentFlux
            (
                surfaceScalarField& phi,
                const volVectorField& U,
                const MRFZones& mrfZones
            ) const;

            //- Reconstruct velocity for transient solvers (after pressure
            //  equation and flux reconstruction).
            void reconstructTransientVelocity
            (
                volVectorField& U,
                surfaceScalarField& phi,
                const fvVectorMatrix& ddtUEqn,
                const volScalarField& rAU,
                const volScalarField& p
            ) const;

            //- Reconstruct velocity for steady state solvers (after pressure
            //  equation and flux reconstruction).
            void reconstructSteadyVelocity
            (
                volVectorField& U,
                const volScalarField& rAU,
                const volScalarField& p
            ) const;

            //- Const access to aCoeff (needed for pressure equation) given the
            //  name of the velocity field
            const surfaceScalarField& aCoeff(const word& UName) const;


        // Evolution

            //- Main control loop
            virtual bool loop() = 0;

            //- Non-orthogonal corrector loop
            inline bool correctNonOrthogonal();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "solutionControlI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "solutionControlTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
