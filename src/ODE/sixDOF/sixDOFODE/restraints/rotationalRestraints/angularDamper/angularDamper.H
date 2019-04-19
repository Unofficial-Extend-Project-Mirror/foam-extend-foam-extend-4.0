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
    Foam::angularDamper

Description
    Rotational restraint corresponding to angular damper defined by angular
    damping coefficients.

Author
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    angularDamper.C

\*---------------------------------------------------------------------------*/

#ifndef angularDamper_H
#define angularDamper_H

#include "rotationalRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class angularDamper Declaration
\*---------------------------------------------------------------------------*/

class angularDamper
:
    public rotationalRestraint
{
    // Private Data

        //- Angular damping coefficients
        diagTensor angDampingCoeffs_;

        //- Whether the damper is applied in global or local c. s.
        Switch inGlobal_;


public:

    //- Runtime type information
    TypeName("angularDamper");


    // Constructors

        //- Construct from dictionary
        angularDamper
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<rotationalRestraint> clone() const;


    // Destructor

        virtual ~angularDamper();


    // Member Functions

        // Restraint specific functions

            //- Return restraining moment (in the global coordinate system)
            virtual vector restrainingMoment
            (
                const scalar t,
                const tensor& toRelative,
                const vector& omega
            ) const;


    // I-O Functions and Operators

        //- Virtual write function
        virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
