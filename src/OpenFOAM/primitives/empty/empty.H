/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    Foam::empty

Description
    A class without storage. Used, for example, in HashSet.

SourceFiles
    emptyI.H
    empty.C
    emptyIO.C

\*---------------------------------------------------------------------------*/

#ifndef empty_H
#define empty_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators

class empty;

Istream& operator>>(Istream& is, empty&);
Ostream& operator<<(Ostream& os, const empty&);


/*---------------------------------------------------------------------------*\
                           Class empty Declaration
\*---------------------------------------------------------------------------*/

class empty
{

public:

    // Constructors

        //- Construct null
        empty()
        {}

        //- Construct from Istream
        empty(Istream&)
        {}


    // IOstream Operator

        friend Istream& operator>>(Istream& is, empty&)
        {
            return is;
        }

        friend Ostream& operator<<(Ostream& os, const empty&)
        {
            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
