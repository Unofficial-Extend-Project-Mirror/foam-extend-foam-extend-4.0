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
    aravasMises

Description
    Stress return method based on Aravas, 1987, On the numerical
    integration of pressure-dependent plasticity models,
    International journal of numerical methods in engineering.
    This is a simplified version for pressure independent plasticity.

SourceFiles
    aravasMises.C

Author
    David Mc Auliffe UCD
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef aravasMises_H
#define aravasMises_H

#include "plasticityStressReturn.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

  class constitutiveModel;

/*---------------------------------------------------------------------------*\
                         Class aravasMises Declaration
\*---------------------------------------------------------------------------*/

class aravasMises
:
    public plasticityStressReturn
{
    // Private data

        // Reference to constitutive model
        constitutiveModel& constitutiveModel_;

        // Yield stress
        volScalarField sigmaY_;

        // Incremental change in sigmaY
        volScalarField DSigmaY_;

        // Increment of plastic strain
        volSymmTensorField DEpsilonP_;

        // Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        // Active yielding flag
        // 1.0 for active yielding
        // 0.0 otherwise
        volScalarField activeYield_;

        // plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        // tolerance and max iterations for Newton loop
        static scalar LoopTol_;
        static label MaxNewtonIter_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        aravasMises(const aravasMises&);

        //- Disallow default bitwise assignment
        void operator=(const aravasMises&);

        // Aravas method functions
        scalar qfun (const scalar qe, const scalar deq, const scalar G) const;

        scalar s0fun (const scalar ebar, const label cellID) const;

        scalar ebarfun
        (
            const scalar ebart,
            const scalar q,
            const scalar deq,
            const scalar s0
        ) const;

        // flow potential
        scalar gfun
    (
     const scalar ebart,
     const scalar deq,
     const scalar qe,
     const scalar G,
     const label cellID
     ) const;

        // Calculate plastic equivalent strain increment and new yield stress
        void aravasNewtonLoop
    (
     scalar& deq,
     scalar& s0,
     const scalar s00,
     const scalar G,
     const scalar ebart,
     const scalar qe,
     const label cellID,
     const scalar maxMagDEpsilon
     ) const;

public:

    //- Runtime type information
    TypeName("aravasMises");

    // Static data members


    // Constructors

        //- Construct from dictionary
        aravasMises
        (
     const word& name,
     constitutiveModel& constitutiveModel
        );


    // Destructor

        virtual ~aravasMises();


    // Member Functions

        //- Correct the rheological model
        virtual void correct();

        //- Return slave friction value fraction
        virtual const volSymmTensorField& DEpsilonP() const
        {
      return DEpsilonP_;
    }

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
