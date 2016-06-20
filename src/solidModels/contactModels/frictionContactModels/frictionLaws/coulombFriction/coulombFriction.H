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
    coulombFriction

Description
    coulombFriction friction law
    shear traction is linearly related to the contact pressure
    via the friction coefficient

SourceFiles
    coulombFriction.C

Author
    Philip Cardiff UCD

\*---------------------------------------------------------------------------*/

#ifndef coulombFriction_H
#define coulombFriction_H

#include "frictionLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class coulombFriction Declaration
\*---------------------------------------------------------------------------*/

class coulombFriction
:
    public frictionLaw
{
    // Private data

         //- friction law coefficientdict
        dictionary frictionLawDict_;

        //- Coulomb friction coefficient
        scalar frictionCoeff_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        coulombFriction(const coulombFriction&);

        //- Disallow default bitwise assignment
        void operator=(const coulombFriction&);

public:

    //- Runtime type information
    TypeName("coulomb");

    // Static data members


    // Constructors

        //- Construct from dictionary
        coulombFriction
        (
            const word& name,
            const dictionary& dict
        );


    // Destructor

        virtual ~coulombFriction();


    // Member Functions

        //- Return slip traction field
        // the traction that would cause slipping
        // for the given pressure
        virtual scalar slipTraction(const scalar pressure);

        //- Write law dictionary
        virtual void writeDict(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
