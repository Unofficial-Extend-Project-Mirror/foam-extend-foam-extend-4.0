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
    geometricSixDOF

Description
    6-DOF solver using a geometric method for integration of rotations.

    Reference (bibtex):

    @article {mullerTerze2016,
        Author = {M\"{u}ller, A. and Terze, Z.},
        title = {Geometric methods and formulations in computational
                 multibody systems},
        Journal = {Acta Mechanica},
        Year = {2016},
        Volume = {227},
        Number = {12},
        Pages = {3327--3350}
    }

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    geometricSixDOF.C

\*---------------------------------------------------------------------------*/

#ifndef geometricSixDOF_H
#define geometricSixDOF_H

#include "sixDOFODE.H"
#include "finiteRotation.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class geometricSixDOF Declaration
\*---------------------------------------------------------------------------*/

class geometricSixDOF
:
    public sixDOFODE
{
    // Private data

        // Initial body state variables

            //- Displacement relative to spring equilibrium
            dimensionedVector Xrel_;

            //- Velocity of mass centroid
            dimensionedVector U_;

            //- Average velocity of mass centroid at previous time-step
            dimensionedVector Uaverage_;

            //- Finite rotation. Note: changed during solution process
            mutable finiteRotation rotation_;

            //- Rotational velocity about mass centroid
            dimensionedVector omega_;


        // Average variables that need to be stored

            //- Average rotational velocity in relative coordinate system
            dimensionedVector omegaAverage_;

            //- Average rotational velocity in absolute coordinate system
            dimensionedVector omegaAverageAbsolute_;


        // ODE controls

            //- Number of equations (depending on rotational constraints)
            scalar nEqns_;

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

            //- Constraints in referent coordinate system
            Switch referentMotionConstraints_;

            //- Rotation of referent coordinate system. Currently defined using
            //  quaternions for comatibility. Need to reformulate.
            //  VV, 28/Feb/2017.
            HamiltonRodriguezRot referentRotation_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        geometricSixDOF(const geometricSixDOF&);

        //- Disallow default bitwise assignment
        void operator=(const geometricSixDOF&);


        // Variables in relative coordinate system (solved for)

            //- Return acceleration in relative coordinate system
            //  given current values of relative displacement and velocity
            dimensionedVector A
            (
                const dimensionedVector& xR,
                const dimensionedVector& uR,
                const tensor& R
            ) const;


            //- Return rotational acceleration in relative coordinate system
            //  given current values for relative rotational velocity
            dimensionedVector OmegaDot
            (
                const tensor& R,
                const dimensionedVector& omega
            ) const;

            //- Return the Euler part of moment equation
            dimensionedVector E
            (
                const dimensionedVector& omega
            ) const;

            //- Constrain translation vector in referent or global coordinate
            //  system
            void constrainTranslation(vector& vec) const;

            //- Exponential map used to calculate increment of the rotation
            //  tensor
            tensor expMap(const vector& rotInc) const;

            //- Differential of the expontential map used to calculate the time
            //  derivative of rotation increment vector
            vector dexpMap(const vector& rotInc, const vector& omega) const;


public:

    // Run-time type information
    TypeName("geometricSixDOF");


    // Static data members

        //- Rotational increment tensor tolerance. Used in expMap member
        //  function in case the rotation is negligibly small
        static const debug::tolerancesSwitch rotIncTensorTol_;

        //- Rotational increment rate of change tolerance. Used in dexpMap
        //  member function in case the rotation rate is negligibly small
        static const debug::tolerancesSwitch rotIncRateTol_;


    // Constructors

        //- Construct from dictionary
        geometricSixDOF(const IOobject& io);

        //- Construct geometricSixDOF object, changing name
        geometricSixDOF
        (
            const word& name,
            const geometricSixDOF& gsd
        );

        //- Return a clone, changing the name
        virtual autoPtr<sixDOFODE> clone(const word& name) const;


    // Destructor

        virtual ~geometricSixDOF();


    // Member Functions

        // Virtual interface for 6DOF motion state

            // Variables in relative coordinate system

                //- Return displacement in translated coordinate system
                //  relative to spring equilibrium
                virtual const dimensionedVector& Xrel() const;

                //- Return rotational velocity in relative coordinate system
                virtual const dimensionedVector& omega() const;


            // Displacement and rotation in the absolute coordinate system

                //- Return position of origin in absolute coordinate system
                virtual dimensionedVector X() const;

                //- Return velocity of origin
                virtual const dimensionedVector& U() const;

                //- Return average velocity of origin for the previous time-step
                virtual const dimensionedVector& Uaverage() const;

                //- Return finite rotation
                virtual const finiteRotation& rotation() const;

                //- Return rotational vector of body
                virtual vector rotVector() const;

                //- Return rotation angle of body
                virtual dimensionedScalar rotAngle() const;


            // Non-access control for setting state variables

                //- Set position of origin
                virtual void setXrel(const vector& x);

                //- Set velocity of origin
                virtual void setU(const vector& u);

                //- Set rotational angle in relative coordinate system
                virtual void setRotation(const HamiltonRodriguezRot& rot);

                //- Set rotational velocity in relative coordinate system
                virtual void setOmega(const vector& omega);

                //- Set referent coordinate system to apply constraints
                virtual void setReferentRotation
                (
                    const HamiltonRodriguezRot& rot
                );

                //- Set ODE parameters from another ODE
                virtual void setState(const sixDOFODE&);


            // Average motion per time-step

                //- Return average rotational vector of body
                virtual vector rotVectorAverage() const;

                //- Return average rotational velocity in relative coordinate
                //  system
                virtual const dimensionedVector& omegaAverage() const;

                //- Return average rotational velocity in absolute coordinate
                //  system
                virtual const dimensionedVector& omegaAverageAbsolute() const;


            // Rotations

                //- Return rotation tensor to relative coordinate system
                virtual tensor toRelative() const;

                //- Return rotation tensor to absolute coordinate system
                virtual tensor toAbsolute() const;

                //- Return transformation tensor between new and previous
                //  rotation
                virtual const tensor& rotIncrementTensor() const;


            // ODE parameters

                //- Return number of equations
                virtual label nEqns() const
                {
                    return nEqns_;
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

                //- Evaluate derivatives
                virtual void derivatives
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dydx
                ) const;

                //- Evaluate Jacobian
                virtual void jacobian
                (
                    const scalar x,
                    const scalarField& y,
                    scalarField& dfdx,
                    scalarSquareMatrix& dfdy
                ) const
                {
                    notImplemented
                    (
                        "geometricSixDOF::jacobian\n"
                        "(\n"
                        "    const scalar x,\n"
                        "    const scalarField& y,\n"
                        "    scalarField& dfdx,\n"
                        "    scalarSquareMatrix& dfdy,\n"
                        ") const"
                    );
                }

                //- Update ODE after the solution, advancing by delta
                virtual void update(const scalar delta);

        // Write controls

            //- WriteData member function required by regIOobject
            virtual bool writeData(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
