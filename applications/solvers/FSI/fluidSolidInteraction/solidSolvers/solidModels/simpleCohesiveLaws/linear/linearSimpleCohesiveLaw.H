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
    linearSimpleCohesiveLaw

Description
    Linear simple cohesive law.

SourceFiles
    linearSimpleCohesiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearSimpleCohesiveLaw_H
#define linearSimpleCohesiveLaw_H

#include "simpleCohesiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearSimpleCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearSimpleCohesiveLaw
:
    public simpleCohesiveLaw
{

private:

    // Private data

        //- Critical separation distance
        dimensionedScalar deltaC_;

public:

    //- Runtime type information
    TypeName("linear");


    // Constructors

        //- Construct from components
        linearSimpleCohesiveLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict
        );

        //- Construct as copy
        linearSimpleCohesiveLaw(const linearSimpleCohesiveLaw&);

        //- Construct and return a clone
        virtual autoPtr<simpleCohesiveLaw> clone() const
        {
            return autoPtr<simpleCohesiveLaw>
            (
                new linearSimpleCohesiveLaw(*this)
            );
        }


    // Destructor

        virtual ~linearSimpleCohesiveLaw();


    // Member Functions

        //- Return reference to critical separation distance
        const dimensionedScalar& deltaC() const
        {
            return deltaC_;
        }

        //- Return current holding traction
        scalar traction(scalar delta) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
