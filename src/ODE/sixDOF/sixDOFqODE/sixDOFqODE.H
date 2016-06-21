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
    sixDOFqODE

Description
    6-DOF solver using quaternions

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Updated by Hrvoje Jasak.

SourceFiles
    sixDOFqODEI.H
    sixDOFqODE.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFqODE_H
#define sixDOFqODE_H

#include "ODE.H"
#include "IOdictionary.H"
#include "dimensionedTypes.H"
#include "HamiltonRodriguezRot.H"
#include "finiteRotation.H"
#include "vectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class sixDOFqODE Declaration
\*---------------------------------------------------------------------------*/

class sixDOFqODE
:
    public IOdictionary,
    public ODE
{
    // Private data

        // Body data

            //- Mass
            dimensionedScalar mass_;

            //- Rotational moment of inertia around centre of mass
            //  in body (relative) coordinates - aligned with main axes
            dimensionedDiagTensor momentOfInertia_;


        // Platform variables

            //- Spring equilibrium position for translation
            dimensionedVector Xequilibrium_;

            //- Linear spring coeffs
            dimensionedDiagTensor linSpringCoeffs_;

            //- Linear damping coeffs
            dimensionedDiagTensor linDampingCoeffs_;

            //- Translational relaxation factor
            scalar relaxFactorT_;

            //- Rotational relaxation factor
            scalar relaxFactorR_;

            //- Old translational relaxation factor
            scalar oldRelaxFactorT_;

            //- Old rotational relaxation factor
            scalar oldRelaxFactorR_;


        // Initial body state variables

            //- Displacement relative to spring equilibrium
            dimensionedVector Xrel_;

            //- Velocity of mass centroid
            dimensionedVector U_;

            //- Average velocity of mass centroid at previous time-step
            dimensionedVector Uaverage_;

            //- Finite rotation
            finiteRotation rotation_;

            //- Rotational velocity about mass centroid
            dimensionedVector omega_;


        // Average variables that need to be stored

            //- Average rotational velocity in relative coordinate system
            dimensionedVector omegaAverage_;

            //- Average rotational velocity in absolute coordinate system
            dimensionedVector omegaAverageAbsolute_;


        // External forces

            //- Accelerations from previous iterations
            //  A_[2] is the old value, A_[1] old old, A_[0] old old old
            List<vector> A_;
            List<vector> OmegaDot_;

            //- Previos iteration non relaxed accelerations
            List<vector> An_;
            List<vector> OmegaDotn_;

            //- Force driving the motion in inertial coord. sys.
            dimensionedVector force_;

            //- Moment driving the motion in inertial coord. sys.
            dimensionedVector moment_;

            //- Force driving the motion in relative coord. sys.
            dimensionedVector forceRelative_;

            //- Moment driving the motion in relative coord. sys.
            dimensionedVector momentRelative_;


        //- ODE coefficients
        scalarField coeffs_;


        //- Motion constraints (given as fixed motion components)

            //- Fixed surge (x-translation)
            Switch fixedSurge_;

            //- Fixed sway (y-translation)
            Switch fixedSway_;

            //- Fixed heave (z-translation)
            Switch fixedHeave_;

            //- Fixed roll (rotation around x)
            Switch fixedRoll_;

            //- Fixed pitch (rotation around y)
            Switch fixedPitch_;

            //- Fixed yaw (rotation around z)
            Switch fixedYaw_;

            //- Restraints in referent coordinate system
            Switch referentMotionConstraints_;

            //- Rotation of referent coordinate system
            HamiltonRodriguezRot referentRotation_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFqODE(const sixDOFqODE&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFqODE&);


        //- Set ODE coefficients from position and rotation
        void setCoeffs();


        // Variables in relative coordinate system (solved for)

            //- Return acceleration in relative coordinate system
            //  given current values of relative displacement and velocity
            dimensionedVector A
            (
                const dimensionedVector& xR,
                const dimensionedVector& uR,
                const HamiltonRodriguezRot& rotation
            ) const;


            //- Return rotational acceleration in relative coordinate system
            //  given current values for relative rotational velocity
            dimensionedVector OmegaDot
            (
                const HamiltonRodriguezRot& rotation,
                const dimensionedVector& omega
            ) const;

            //- Return the Euler part of moment equation
            dimensionedVector E
            (
                const dimensionedVector& omega
            ) const;

            //- Constrain rotation vector in referent or global coordinate
            //  system
            void constrainRotation(vector& vec) const;

            //- Constrain translation vector in referent or global coordinate
            //  system
            void constrainTranslation(vector& vec) const;

            //- Update Aitkens relaxation parameters
            void aitkensRelaxation(const scalar min, const scalar max);


public:

    // Constructors

        //- Construct from dictionary
        sixDOFqODE(const IOobject& io);

        //- Construct sixDOFqODE, changing name
        sixDOFqODE
        (
            const word& name,
            const sixDOFqODE& sd
        );


    // Destructor

        virtual ~sixDOFqODE();


    // Member Functions

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


        // Variables in relative coordinate system (solved for)

            //- Return displacement in translated coordinate system
            //  relative to spring equilibrium
            inline const dimensionedVector& Xrel() const;

            //- Return rotational velocity in relative coordinate system
            inline const dimensionedVector& omega() const;

            //- Return rotational acceleration in relative coordinate system
            inline dimensionedVector omegaDot() const;


        // Displacement and rotation in the absolute coordinate system

            //- Return position of origin in absolute coordinate system
            inline dimensionedVector X() const;

            //- Return velocity of origin
            inline const dimensionedVector& U() const;

            //- Return average velocity of origin for the previous time-step
            inline const dimensionedVector& Uaverage() const;

            //- Return finite rotation
            inline const finiteRotation& rotation() const;

            //- Return acceleration of origin
            inline dimensionedVector A() const;

            //- Return rotational vector of body
            inline vector rotVector() const;

            //- Return rotation angle of body
            inline dimensionedScalar rotAngle() const;


        // Non-constant access

            //- Set position of origin
            inline void setXrel(const vector& x);

            //- Set velocity of origin
            inline void setU(const vector& u);

            //- Set rotational angle in relative
            //  coordinate system
            inline void setRotation(const HamiltonRodriguezRot& rot);

            //- Set rotational velocity in relative
            //  coordinate system
            inline void setOmega(const vector& omega);

            //- Set referent coordinate system to apply constraints
            inline void setReferentRotation(const HamiltonRodriguezRot& rot);


        // Average motion per time-step

            //- Return average rotational vector of body
            inline vector rotVectorAverage() const;

            //- Return average rotational velocity in
            //  relative coordinate system
            inline const dimensionedVector& omegaAverage() const;

            //- Return average rotational velocity in
            //  absolute coordinate system
            inline const dimensionedVector& omegaAverageAbsolute() const;


        // Force and moment

            //- Return force in inertial coord. sys.
            inline const dimensionedVector& force() const;

            //- Return access to force in inertial coord. sys.
            inline dimensionedVector& force();

            //- Return moment in inertial coord. sys.
            inline const dimensionedVector& moment() const;

            //- Return access to moment in inertial coord. sys.
            inline dimensionedVector& moment();

            //- Return force in relative coord. sys.
            inline const dimensionedVector& forceRelative() const;

            //- Return access to force in relative coord. sys.
            inline dimensionedVector& forceRelative();

            //- Return moment in relative coord. sys.
            inline const dimensionedVector& momentRelative() const;

            //- Return access to moment in relative coord. sys.
            inline dimensionedVector& momentRelative();

            //- Return total force in inertial coord. sys.
            inline dimensionedVector forceTotal() const;

            //- Return total moment in inertial coord. sys.
            inline dimensionedVector momentTotal() const;

            //- Relax the force(acceleration) using Aitkens or fixed relaxation
            void relaxAcceleration
            (
                const scalar minRelFactor,
                const scalar maxRelFactor
            );


        // Rotations

            //- Return rotation tensor to relative coordinate system
            inline tensor toRelative() const;

            //- Return rotation tensor to absolute coordinate system
            inline tensor toAbsolute() const;

            //- Return transformation tensor between new and previous rotation
            inline const tensor& rotIncrementTensor() const;

            //-  Set ODE parameters from another ODE
            void setState(const sixDOFqODE&);


        // ODE parameters

            //- Return number of equations
            virtual label nEqns() const
            {
                return 13;
            }

            //- Return access to coefficients
            virtual scalarField& coeffs()
            {
                return coeffs_;
            }

            //- Return reference to coefficients
            virtual const scalarField& coeffs() const
            {
                return coeffs_;
            }

            //- Return derivatives
            virtual void derivatives
            (
                const scalar x,
                const scalarField& y,
                scalarField& dydx
            ) const;

            //- Return Jacobian
            virtual void jacobian
            (
                const scalar x,
                const scalarField& y,
                scalarField& dfdx,
                scalarSquareMatrix& dfdy
            ) const;

            //- Update ODE after the solution, advancing by delta
            virtual void update(const scalar delta);


        //- WriteData member function required by regIOobject
        bool writeData(Ostream&) const;


    // Ostream operator

        friend Ostream& operator<<(Ostream&, const sixDOFqODE&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "sixDOFqODEI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
