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
    Foam::labelledTri

Description
    Triangle with additional region number.

SourceFiles
    labelledTriI.H

\*---------------------------------------------------------------------------*/

#ifndef labelledTri_H
#define labelledTri_H

#include "triFace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class labelledTri Declaration
\*---------------------------------------------------------------------------*/

class labelledTri
:
    public triFace
{
    // Private data

        //- Region ID
        label region_;


public:

    // Constructors

        //- Construct null
        inline labelledTri();

        //- Construct from three point labels and a region label
        inline labelledTri
        (
            const label a,
            const label b,
            const label c,
            const label region
        );

        //- Construct from Istream
        inline labelledTri(Istream&);


    // Member Functions

        // Access

            //- Return region label
            inline label region() const;

            //- Return region label
            inline label& region();


        // Check

        // Edit

        // Write


    // Friend Functions

    // Friend Operators

    // IOstream Operators

        inline friend Istream& operator>>(Istream&, labelledTri&);
        inline friend Ostream& operator<<(Ostream&, const labelledTri&);
};


template<>
inline bool contiguous<labelledTri>()  {return true;}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "labelledTriI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
