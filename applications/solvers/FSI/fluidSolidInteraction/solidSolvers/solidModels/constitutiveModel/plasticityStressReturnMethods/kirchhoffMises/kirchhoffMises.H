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
    kirchhoffMises

Description
    J2 with isotropic hardening return mapping for Kirchhoff stress based
    Updated Lagrangian elasto-plastic solver.

    Procedure similar to approach outlined in Box 9.1 in Simo & Hughes,
    Computational Inelasticity, 1998, Springer.

    If plasticity is nonlinear Newton's method is used to calculate the
    plastic strain increment.

SourceFiles
    kirchhoffMises.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef kirchhoffMises_H
#define kirchhoffMises_H

#include "plasticityStressReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class constitutiveModel;

/*---------------------------------------------------------------------------*\
                         Class kirchhoffMises Declaration
\*---------------------------------------------------------------------------*/

class kirchhoffMises
:
    public plasticityStressReturn
{
    // Private data

        // Reference to constitutive model
        constitutiveModel& constitutiveModel_;

        // Restart case
        bool restarted_;

        // Yield stress
        volScalarField sigmaY_;
        surfaceScalarField sigmaYf_;

        // Incremental change in sigmaY
        volScalarField DSigmaY_;
        surfaceScalarField DSigmaYf_;

        // Increment of plastic strain volField
        volSymmTensorField DEpsilonP_;
        surfaceSymmTensorField DEpsilonPf_;

//         // Elastic left Cauchy-Green strain tensor
//         // with volumetric term removed
//         volSymmTensorField bEbar_;
//         surfaceSymmTensorField bEbarf_;

        // Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;
        surfaceScalarField DEpsilonPEqf_;

        // Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;
        surfaceScalarField DLambdaf_;

        // Equivalent plastic strain
        volScalarField epsilonPEq_;
        surfaceScalarField epsilonPEqf_;

        // Active yielding flag
        // 1.0 for active yielding
        // 0.0 otherwise
        volScalarField activeYield_;

        // plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;
        surfaceSymmTensorField plasticNf_;

        // An iterative procedure is used when the yield stress is a nonlinear
        // function of plastic strain
        const bool nonLinearPlasticity_;

        // Linear plastic modulus. It is not defined if plastic is nonlinear
        // or if plastic is perfectly plastic
        const volScalarField* HPtr_;
        const surfaceScalarField* HfPtr_;

        // Whether or not to extrapolate DEpsilonP to the boundary, instead of
        // calculating it at the patch. This helps fix convergence problem with
        // boundaries such as contact.
        boolList extrapBoundary_;

        // Tolerance for Newton loop
        static scalar LoopTol_;

        // Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        // finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        // Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        kirchhoffMises(const kirchhoffMises&);

        //- Disallow default bitwise assignment
        void operator=(const kirchhoffMises&);

        // Evaulate current yield stress
        scalar curYieldStress
        (
            const scalar curEpsilonPEq,
            const scalar J, // Current Jacobian
            const label cellID
        ) const;

        // Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar magSTrial,
            const scalar DLambda,
            const scalar muBar,
            const scalar J, // Current Jacobian
            const label cellID
        ) const;

        // Iteratively calculate plastic multiplier increment (DLambda)
        // and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,
            scalar& curSigmaY,
            const scalar magSTrial,
            const scalar muBar,
            const scalar J, // Current Jacobian
            const label cellID,
            const scalar maxMagDEpsilon
        ) const;

        // Iteratively calculate plastic multiplier increment (DLambda)
        // and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,
            scalar& curSigmaY,
            const scalar magSTrial,
            const scalar muBar,
            const scalar epsilonPEq,
            const scalar J, // Current Jacobian
            const label cellID
        ) const;

public:

    //- Runtime type information
    TypeName("kirchhoffMises");

    // Static data members


    // Constructors

        //- Construct from dictionary
        kirchhoffMises
        (
            const word& name,
            constitutiveModel& constitutiveModel
        );


    // Destructor

        virtual ~kirchhoffMises();


    // Member Functions

        //- Correct DEpsilonPf or DEpsilonP
        virtual void correct();

        //- Correct DEpsilonPf
        virtual void correctf();

        //- Correct DEpsilonP
        virtual void correctCell();

        //- Correct DEpsilonPf
        virtual void correctFace();

        //- Return const reference to plastic strain increment
        virtual const volSymmTensorField& DEpsilonP() const;

        //- Return const reference to plastic strain increment
        virtual const surfaceSymmTensorField& DEpsilonPf() const;

        //- Return slave friction value fraction
        virtual void updateYieldStress();

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const {};
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
