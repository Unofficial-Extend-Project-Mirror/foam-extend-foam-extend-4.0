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
    Foam::periodicOscillation

Description
    Periodic translational motion given by:

    g(vDot, t) = vDot + a*omega^2*sin(omega*t + phi) = 0,

    where a is the amplitude of oscillation, omega radial frequency and phi is
    the phase shift.

    Note: since this constraint is basically algebraic (dependent on
    displacement), the differentiation index is two so the stabilisation on
    displacement and velocity is necessary.

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    periodicOscillation.C

\*---------------------------------------------------------------------------*/

#ifndef periodicOscillation_H
#define periodicOscillation_H

#include "translationalConstraint.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class periodicOscillation Declaration
\*---------------------------------------------------------------------------*/

class periodicOscillation
:
    public translationalConstraint
{
    // Private data

        //- Direction
        vector dir_;

        //- Amplitude of the motion
        const scalar a_;

        //- Period of the motion
        const scalar period_;

        //- Radian frequency of the motion
        const scalar omega_;

        //- Phase shift (in degrees)
        const scalar phi_;


public:

    //- Runtime type information
    TypeName("periodicOscillation");


    // Constructors

        //- Construct from dictionary
        periodicOscillation
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<translationalConstraint> clone() const;


    // Destructor

        virtual ~periodicOscillation();


    // Member Functions

        // Constraint specific functions

            //- Return matrix contribution defined by constraint, f(t)
            virtual vector matrixContribution
            (
                const scalar,
                const tensor&,
                const vector&,
                const vector&
            ) const;

            //- Return source contribution defined by constraint, a(t)
            virtual scalar sourceContribution
            (
                const scalar t,
                const tensor&,
                const vector&,
                const vector&
            ) const;

            //- Stabilise the constraint
            virtual void stabilise
            (
                const scalar t,
                vector& x,
                vector& u
            ) const;


    // I-O Functions

        //- Virtual write function
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
