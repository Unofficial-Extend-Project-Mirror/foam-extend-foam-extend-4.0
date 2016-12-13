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
    Foam::STLpoint

Description
    A vertex point representation for STL files.

SourceFiles
    STLpointI.H

\*---------------------------------------------------------------------------*/

#ifndef STLpoint_H
#define STLpoint_H

#include "point.H"
#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class STLpoint Declaration
\*---------------------------------------------------------------------------*/

class STLpoint
:
    public Vector<float>
{

public:

    // Constructors

        //- Construct null
        inline STLpoint()
        {}

        //- Construct from components
        inline STLpoint(float x, float y, float z)
        :
            Vector<float>(x, y, z)
        {}

        //- Construct from point
        inline STLpoint(const point& pt)
        :
            Vector<float>(float(pt.x()), float(pt.y()), float(pt.z()))
        {}

        //- Construct from istream
        inline STLpoint(Istream& is)
        :
            Vector<float>(is)
        {}


    // Member Operators

        //- Conversion to point
        inline operator point() const
        {
            return point(x(), y(), z());
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
