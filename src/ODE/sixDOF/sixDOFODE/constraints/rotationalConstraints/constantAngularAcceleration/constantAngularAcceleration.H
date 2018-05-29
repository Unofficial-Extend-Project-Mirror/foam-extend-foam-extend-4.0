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
    Foam::constantAngularAcceleration

Description
    Rotational constraint defined by constant angular acceleration:

    g(omegaDot) = omegaDot - a = 0.

Author
    Viktor Pandza, FSB Zagreb.  All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    constantAngularAcceleration.C

\*---------------------------------------------------------------------------*/

#ifndef constantAngularAcceleration_H
#define constantAngularAcceleration_H

#include "rotationalConstraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class constantAngularAcceleration Declaration
\*---------------------------------------------------------------------------*/

class constantAngularAcceleration
:
    public rotationalConstraint
{
    // Private data

        //- Direction of the constraint (unit vector)
        vector dir_;

        //- Constant value of angular acceleration
        const scalar alpha_;

        //- Switch whether the constraint should be applied in local or global
        //  coordinate system
        const Switch inGlobal_;


public:

    //- Runtime type information
    TypeName("constantAngularAcceleration");


    // Constructors

        //- Construct from dictionary
        constantAngularAcceleration
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<rotationalConstraint> clone() const;


    // Destructor

        virtual ~constantAngularAcceleration();


    // Member Functions

        // Constraint specific functions

            //- Return matrix contribution defined by constraint, f(t)
            virtual vector matrixContribution
            (
                const scalar,
                const tensor& toRelative,
                const vector&
            ) const;

            //- Return source contribution defined by constraint, a(t)
            virtual scalar sourceContribution
            (
                const scalar,
                const tensor&,
                const vector&
            ) const;

            //- Stabilise the constraint
            virtual void stabilise
            (
                const scalar,
                vector&
            ) const
            {
                // Does nothing: no need to stabilise this constraint
            }


    // I-O Functions

        //- Virtual write function
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
