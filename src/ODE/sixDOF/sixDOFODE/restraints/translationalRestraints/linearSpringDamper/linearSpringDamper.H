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
    Foam::linearSpringDamper

Description
    Translational restraint corresponding to linear spring and linear damper
    defined by linear spring coefficients and linear damping coefficients.

Author
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    linearSpringDamper.C

\*---------------------------------------------------------------------------*/

#ifndef linearSpringDamper_H
#define linearSpringDamper_H

#include "translationalRestraint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class linearSpringDamper Declaration
\*---------------------------------------------------------------------------*/

class linearSpringDamper
:
    public translationalRestraint
{
    // Private Data

        //- Linear spring coefficients
        diagTensor linSpringCoeffs_;

        //- Linear damping coefficients
        diagTensor linDampingCoeffs_;


public:

    //- Runtime type information
    TypeName("linearSpringDamper");


    // Constructors

        //- Construct from dictionary
        linearSpringDamper
        (
            const word& name,
            const dictionary& dict,
            const sixDOFODE& sixDOF
        );

        //- Construct and return a clone
        virtual autoPtr<translationalRestraint> clone() const;


    // Destructor

        virtual ~linearSpringDamper();


    // Member Functions

        // Restraint specific functions

            //- Return restraining force (in the global coordinate system)
            virtual vector restrainingForce
            (
                const scalar t,
                const tensor& toRelative,
                const vector& x,
                const vector& u
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
