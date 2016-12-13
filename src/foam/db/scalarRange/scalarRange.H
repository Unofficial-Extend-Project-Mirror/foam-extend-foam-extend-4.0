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
    Foam::scalarRange

Description
    A scalar range specifier.

    The range selector can be specified as an "LOWER:UPPER" range, as a
    "LOWER:" bound, as an ":UPPER" bound or simply as an "EXACT" value.
    The read constructor uses a colon (:) as a range marker and a comma (,)
    to delimit the next possible range selector.

SourceFiles
    scalarRange.C

\*---------------------------------------------------------------------------*/
#ifndef scalarRange_H
#define scalarRange_H

#include "scalar.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class Istream;
class Ostream;

// Forward declaration of friend functions and operators
class scalarRange;
Istream& operator>>(Istream&, scalarRange&);
Ostream& operator<<(Ostream&, const scalarRange&);


/*---------------------------------------------------------------------------*\
                        Class scalarRange Declaration
\*---------------------------------------------------------------------------*/

class scalarRange
{
    //- Enumeration defining the types of token
    enum rangeType
    {
        EMPTY = 0,
        EXACT,
        LOWER,
        UPPER,
        RANGE
    };


    // Private data

        enum rangeType type_;
        scalar value_;
        scalar value2_;


public:

        static debug::debugSwitch debug;


    // Constructors

        //- Construct Null
        scalarRange();

        //- Construct a Range
        scalarRange(const scalar& lower, const scalar& upper);

        //- Construct from Istream.
        //  Since commas can be used as list delimiters,
        //  leading and trailing commas are ignored.
        scalarRange(Istream&);


    // Member Functions

        //- Is the range non-empty?
        bool isDefined() const;

        //- Is the range 'EXACT'?
        bool isExact() const;

        //- The value constituting an 'EXACT' match
        //  or the values for 'UPPER' or 'LOWER' limits
        scalar value() const;

        //- The lower limit
        scalar lower() const;

        //- The upper limit
        scalar upper() const;

        //- Return true if the value is within the range
        bool selected(const scalar&) const;


    // Member Operators

        bool operator==(const scalarRange&) const;
        bool operator!=(const scalarRange&) const;


    // IOstream Operators

        friend Istream& operator>>(Istream&, scalarRange&);
        friend Ostream& operator<<(Ostream&, const scalarRange&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
