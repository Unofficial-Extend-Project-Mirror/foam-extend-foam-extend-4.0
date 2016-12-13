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
    Foam::tetherPotentials::restrainedHarmonicSpring

Description


SourceFiles
    restrainedHarmonicSpring.C

\*---------------------------------------------------------------------------*/

#ifndef restrainedHarmonicSpring_H
#define restrainedHarmonicSpring_H

#include "tetherPotential.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace tetherPotentials
{

/*---------------------------------------------------------------------------*\
                  Class restrainedHarmonicSpring Declaration
\*---------------------------------------------------------------------------*/

class restrainedHarmonicSpring
:
    public tetherPotential
{
    // Private data

        dictionary restrainedHarmonicSpringCoeffs_;

        scalar springConstant_;

        scalar rR_;


public:

    //- Runtime type information
    TypeName("restrainedHarmonicSpring");


    // Constructors

        //- Construct from components
        restrainedHarmonicSpring
        (
            const word& name,
            const dictionary& tetherPotentialProperties
        );


    // Destructor

        ~restrainedHarmonicSpring()
        {}


    // Member Functions

        scalar energy(const vector r) const;

        vector force(const vector r) const;

        //- Read dictionary
        bool read(const dictionary& tetherPotentialProperties);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace tetherPotentials
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
