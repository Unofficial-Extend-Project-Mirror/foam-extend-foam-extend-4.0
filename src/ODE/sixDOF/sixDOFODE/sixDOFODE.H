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
    equations with a fat interface in order to provide representation of
    rotations in quaternions or rotation tensors.

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
#include "sixDOFODEIO.H"
#include "dimensionedTypes.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declarations of classes
class IOobject;

class finiteRotation;
class HamiltonRodriguezRot;

/*---------------------------------------------------------------------------*\
                         Class sixDOFODE Declaration
\*---------------------------------------------------------------------------*/

class sixDOFODE
:
    public ODE
{
    // Private data

        //- Dictionary object controlling I/O for sixDOFODE
        const sixDOFODEIO dict_;


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


        // External forces

            //- Force driving the motion in inertial coord. sys.
            dimensionedVector force_;

            //- Moment driving the motion in inertial coord. sys.
            dimensionedVector moment_;

            //- Force driving the motion in relative coord. sys.
            dimensionedVector forceRelative_;

            //- Moment driving the motion in relative coord. sys.
            dimensionedVector momentRelative_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFODE(const sixDOFODE&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFODE&);


protected:

    // Protected Member Functions

        //- Set ODE coefficients from position and rotation
        virtual void setCoeffs() = 0;

        //- Update Aitkens relaxation parameters
        void aitkensRelaxation(const scalar min, const scalar max);


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

        //- Return a clone, changing the name
        virtual autoPtr<sixDOFODE> clone(const word& name) const = 0;


    // Selectors

        //- Return autoPtr to the selected sixDOFODE
        static autoPtr<sixDOFODE> New(const IOobject& io);

    // Destructor

        virtual ~sixDOFODE();


    // Member Functions

        // Access to common data

            //- Return sixDOFODEIO (IOdictionary)
            const sixDOFODEIO& dict() const;

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

            //- Relax the force (acceleration) using Aitkens or fixed relaxation
            void relaxAcceleration
            (
                const scalar minRelFactor,
                const scalar maxRelFactor
            );


        // Virtual interface for 6DOF motion state

            // Variables in relative coordinate system

                //- Return displacement in translated coordinate system
                //  relative to spring equilibrium
                virtual const dimensionedVector& Xrel() const = 0;

                //- Return rotational velocity in relative coordinate system
                virtual const dimensionedVector& omega() const = 0;


            // Displacement and rotation in the absolute coordinate system

                //- Return position of origin in absolute coordinate system
                virtual dimensionedVector X() const = 0;

                //- Return velocity of origin
                virtual const dimensionedVector& U() const = 0;

                //- Return average velocity of origin for the previous time-step
                virtual const dimensionedVector& Uaverage() const = 0;

                //- Return finite rotation
                virtual const finiteRotation& rotation() const = 0;

                //- Return rotational vector of body
                virtual vector rotVector() const = 0;

                //- Return rotation angle of body
                virtual dimensionedScalar rotAngle() const = 0;


            // Non-access control for setting state variables

                //- Set position of origin
                virtual void setXrel(const vector& x) = 0;

                //- Set velocity of origin
                virtual void setU(const vector& u) = 0;

                //- Set rotational angle in relative coordinate system
                virtual void setRotation(const HamiltonRodriguezRot& rot) = 0;

                //- Set rotational velocity in relative coordinate system
                virtual void setOmega(const vector& omega) = 0;

                //- Set referent coordinate system to apply constraints
                virtual void setReferentRotation
                (
                    const HamiltonRodriguezRot& rot
                ) = 0;

                //- Set ODE parameters from another ODE
                virtual void setState(const sixDOFODE&);


            // Average motion per time-step

                //- Return average rotational vector of body
                virtual vector rotVectorAverage() const = 0;

                //- Return average rotational velocity in relative coordinate
                //  system
                virtual const dimensionedVector& omegaAverage() const = 0;

                //- Return average rotational velocity in absolute coordinate
                //  system
                virtual const dimensionedVector&
                omegaAverageAbsolute() const = 0;


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
