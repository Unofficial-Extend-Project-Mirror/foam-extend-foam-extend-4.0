/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    linearCohesiveLaw

Description
    Linear cohesive law.

SourceFiles
    linearCohesiveLaw.C

\*---------------------------------------------------------------------------*/

#ifndef linearCohesiveLaw_H
#define linearCohesiveLaw_H

#include "cohesiveLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class linearCohesiveLaw Declaration
\*---------------------------------------------------------------------------*/

class linearCohesiveLaw
:
    public cohesiveLaw
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
        linearCohesiveLaw
        (
            const word& cohesiveLawName,
            const dictionary& dict
        );

        //- Construct as copy
        linearCohesiveLaw(const linearCohesiveLaw&);

        //- Construct and return a clone
        virtual autoPtr<cohesiveLaw> clone() const
        {
            return autoPtr<cohesiveLaw>(new linearCohesiveLaw(*this));
        }


    // Destructor

        virtual ~linearCohesiveLaw();


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
