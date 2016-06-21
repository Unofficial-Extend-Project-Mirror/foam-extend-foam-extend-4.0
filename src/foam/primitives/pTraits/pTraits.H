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
    Foam::pTraits

Description
    Traits class for primitives.

    All primitives need a specialised version of this class

\*---------------------------------------------------------------------------*/

#ifndef pTraits_H
#define pTraits_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;

/*---------------------------------------------------------------------------*\
                           Class pTraits Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class pTraits
:
    public PrimitiveType
{

public:

    // Constructors

        pTraits(const PrimitiveType& p)
        :
            PrimitiveType(p)
        {}

        pTraits(Istream& is)
        :
            PrimitiveType(is)
        {}


    // Member operators

    //HJ, unnecessary compilation problems.  HJ, 28/Sep/2010
//         operator PrimitiveType() const
//         {
//             return *this;
//         }

//         operator PrimitiveType&()
//         {
//             return *this;
//         }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
