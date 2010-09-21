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
            const dimensionedVector Xequilibrium_;

            //- Linear spring coeffs
            const dimensionedDiagTensor linSpringCoeffs_;

            //- Linear damping coeffs
            const dimensionedDiagTensor linDampingCoeffs_;


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


public:


//     //- Runtime type information
//       Not possible because of I/O error: incorrect type, expecting dictionary
//       HJ, 11/Feb/2008
//     TypeName("sixDOFqODE");


    // Constructors

        //- Construct from dictionary
        sixDOFqODE(const IOobject& io);


    // Destructor

        virtual ~sixDOFqODE();


    // Member Functions

        //- Return mass
        inline const dimensionedScalar& mass() const;

        //- Return moment of inertia
        inline const dimensionedDiagTensor& momentOfInertia() const;

        //- Return equilibrium position of origin
        inline const dimensionedVector& Xequilibrium() const;


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

            //- Return acceleration of origin
            inline dimensionedVector A() const;

            //- Return rotational vector of body
            inline vector rotVector() const;

            //- Return rotation angle of body
            inline dimensionedScalar rotAngle() const;


        // Non-constant access to velocities

            //- Return access to velocity of origin
            inline dimensionedVector& U();

            //- Return access to rotational velocity in relative
            //  coordinate system
            inline dimensionedVector& omega();


        // Average motion per time-step

            //- Return average rotational vector of body
            inline vector rotVectorAverage() const;

            //- Return average rotational velocity in relative coordinate system
            inline const dimensionedVector& omegaAverage() const;

            //- Return average rotational velocity in absolute coordinate system
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


        // Rotations

            //- Return rotation tensor to relative coordinate system
            inline tensor toRelative() const;

            //- Return rotation tensor to absolute coordinate system
            inline tensor toAbsolute() const;

            //- Return transformation tensor between new and previous rotation
            inline const tensor& rotIncrementTensor() const;


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
