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
    sixDOFODE

Description
    Abstract base class for six-degrees-of-freedom (6DOF) ordinary differential
    equations with necessary interface for two-way coupling with CFD solvers.

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    sixDOFODEI.H
    sixDOFODE.C
    newSixDOFODE.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFODE_H
#define sixDOFODE_H

#include "ODE.H"
#include "dimensionedTypes.H"
#include "autoPtr.H"
#include "OutputControlDictionary.H"
#include "runTimeSelectionTables.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sixDOFODE Declaration
\*---------------------------------------------------------------------------*/

class sixDOFODE
:
    public ODE
{
    // Private data

        //- Dictionary object controlling I/O for sixDOFODE
        const OutputControlDictionary<sixDOFODE> dict_;


        // Body data

            //- Mass
            dimensionedScalar mass_;

            //- Rotational moment of inertia around centre of mass
            //  in body (relative) coordinates - aligned with main axes
            dimensionedDiagTensor momentOfInertia_;


        // Equilibrium position and spring/damping coefficients

            //- Spring equilibrium position for translation
            dimensionedVector Xequilibrium_;

            //- Linear spring coeffs
            dimensionedDiagTensor linSpringCoeffs_;

            //- Linear damping coeffs
            dimensionedDiagTensor linDampingCoeffs_;


        // Aitkens relaxation data members

            //- Switch to control whether to use Aitkens relaxation
            Switch aitkensRelaxation_;

            //- Minimum acceleration relaxation factor (default 0.1)
            const scalar minRelaxFactor_;

            //- Maximum acceleration relaxation factor (default 0.5)
            const scalar maxRelaxFactor_;

            //- Translational relaxation factor
            scalar relaxFactorT_;

            //- Rotational relaxation factor
            scalar relaxFactorR_;

            //- Old translational relaxation factor
            scalar oldRelaxFactorT_;

            //- Old rotational relaxation factor
            scalar oldRelaxFactorR_;

            //- Accelerations from previous iterations
            //  A_[2] is the old value, A_[1] old old, A_[0] old old old
            List<vector> A_;
            List<vector> OmegaDot_;

            //- Previos iteration non relaxed accelerations
            List<vector> An_;
            List<vector> OmegaDotn_;


        // External forces and moments

            //- Force driving the motion in global (inertial) coord. sys.
            dimensionedVector force_;

            //- Moment driving the motion in global (inertial) coord. sys.
            dimensionedVector moment_;


        // Private data used to control multiple ODE updates per time step
        // Note: Before solving the ODE from the top level, we will store the
        // previous state if this is the first update in a given time step. The
        // state will be kept as the pointer to sixDOFODE, which is not optimal
        // because we are keeping track of some unnecessary data. We could use
        // more elegant and efficient code design. VV, 1/Mar/2017.

            //- Local time index
            label curTimeIndex_;

            //- Pointer to the sixDOFODE object carrying old state
            autoPtr<sixDOFODE> oldStatePtr_;


    // Private Member Functions

        // Copy control

            //- Disallow default bitwise copy construct
            sixDOFODE(const sixDOFODE&);

            //- Disallow default bitwise assignment
            void operator=(const sixDOFODE&);


        // Aitkens relaxation helper function

            //- Update Aitkens relaxation parameters
            void updateRelaxFactors();

            //- Relax the force (acceleration) using Aitkens or fixed relaxation
            void relaxAcceleration();


        // Helper function for controlling multiple ODE solver calls per
        // time-step

            //- Initialise ODE before setting external forces/moments and
            //  solving
            virtual void initState();


protected:

    // Protected Member Functions

        // Non-access control for setting state variables

            //- Set ODE parameters from another ODE
            virtual void setState(const sixDOFODE&);


public:

    //- Run-time type information
    TypeName("sixDOFODE");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            sixDOFODE,
            dictionary,
            (const IOobject& io),
            (io)
        );


    // Constructors

        //- Construct from dictionary
        sixDOFODE(const IOobject& io);

        //- Copy construct given new name
        sixDOFODE(const word& name, const sixDOFODE& sd);

        //- Return a clone, changing the name
        virtual autoPtr<sixDOFODE> clone(const word& name) const = 0;


    // Selectors

        //- Return autoPtr to the selected sixDOFODE
        static autoPtr<sixDOFODE> New(const IOobject& io);


    // Destructor

        virtual ~sixDOFODE();


    // Member Functions

        // Access to common data

            //- Return write controlled dictionary
            const OutputControlDictionary<sixDOFODE>& dict() const;

            //- Return mass
            inline const dimensionedScalar& mass() const;

            //- Return access to mass
            inline dimensionedScalar& mass();

            //- Return moment of inertia
            inline const dimensionedDiagTensor& momentOfInertia() const;

            //- Return access to moment of inertia
            inline dimensionedDiagTensor& momentOfInertia();

            //- Return equilibrium position of origin
            inline const dimensionedVector& Xequilibrium() const;

            //- Return access to equilibrium position of origin
            inline dimensionedVector& Xequilibrium();

            //- Return linear spring coefficient
            inline const dimensionedDiagTensor& linSpringCoeffs() const;

            //- Return linear damping coefficient
            inline const dimensionedDiagTensor& linDampingCoeffs() const;


        // Access to forces and moments

            //- Return force in global (inertial) coord. sys.
            inline const dimensionedVector& force() const;

            //- Return moment in global (inertial) coord. sys.
            inline const dimensionedVector& moment() const;

            //- Set external force and moment
            inline void setExternalForceAndMoment
            (
                const vector& externalForce,
                const vector& externalMoment
            );


        // Virtual interface for 6DOF motion state

            // Variables in relative coordinate system

                //- Return displacement in translated coordinate system
                //  relative to spring equilibrium
                virtual const dimensionedVector& Xrel() const = 0;

                //- Return rotational velocity in relative coordinate system
                virtual const dimensionedVector& omega() const = 0;


            // Displacement and velocity in the absolute coordinate system

                //- Return position of origin in absolute coordinate system
                virtual dimensionedVector X() const = 0;

                //- Return velocity of origin
                virtual const dimensionedVector& U() const = 0;

                //- Return average velocity of origin (evaluated at midstep)
                virtual const dimensionedVector& Uaverage() const = 0;


            // Average motion per time-step

                //- Return average rotational velocity in relative coordinate
                //  system (evaluated at midstep)
                virtual const dimensionedVector& omegaAverage() const = 0;


            // Rotations

                //- Return rotation tensor to relative coordinate system
                virtual tensor toRelative() const = 0;

                //- Return rotation tensor to absolute coordinate system
                virtual tensor toAbsolute() const = 0;

                //- Return transformation tensor between new and previous
                //  rotation
                virtual const tensor& rotIncrementTensor() const = 0;


            // ODE parameters

                //- Return number of equations
                virtual label nEqns() const = 0;

                //- Return access to coefficients
                virtual scalarField& coeffs() = 0;

                //- Return reference to coefficients
                virtual const scalarField& coeffs() const = 0;

                //- Evaluate derivatives
                virtual void derivatives
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dydx
                ) const = 0;

                //- Evaluate Jacobian
                virtual void jacobian
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dfdx,
                    scalarSquareMatrix& dfdy
                ) const = 0;

                //- Update ODE after the solution, advancing by delta
                virtual void update(const scalar delta) = 0;


        // Write control

            //- writeData member function required by regIOobject
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDOFODEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
