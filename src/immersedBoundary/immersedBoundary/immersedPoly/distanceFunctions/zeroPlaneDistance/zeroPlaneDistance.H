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
    Foam::zeroPlaneDistance

Description
    Distance class to zero plane

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef zeroPlaneDistance_H
#define zeroPlaneDistance_H

#include "scalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class zeroPlaneDistance Declaration
\*---------------------------------------------------------------------------*/

class zeroPlaneDistance
{
public:

    // Constructors

        //- Construct null
        zeroPlaneDistance();


    //- Destructor - default


    // Member Functions

        // Iterate distance
        inline bool iterateDistance() const
        {
            return false;
        }

        // Return distance to xy-plane at z = 0
        inline scalar distance(const point& p) const
        {
            return p.z();
        }

        // Return distance to xy-plane at z = 0
        inline tmp<scalarField> distance(const vectorField& p) const
        {
            return p.component(vector::Z);
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
