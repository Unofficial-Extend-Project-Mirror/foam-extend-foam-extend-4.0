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
    Foam::safeBool

Description
    Implementation of the Safe Bool Idiom.
    This is the version without virtual functions.

    Author: Bjorn Karlsson
    http://www.artima.com/cppsource/safebool.html

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef safeBool_H
#define safeBool_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class safeBool_base Declaration
\*---------------------------------------------------------------------------*/

class safeBoolBase
{
    public:
        typedef void (safeBoolBase::*boolType)() const;
        void thisTypeDoesNotSupportComparisons() const {}

    protected:

        safeBoolBase() {}
        safeBoolBase(const safeBoolBase&) {}
        safeBoolBase& operator=(const safeBoolBase&) {return *this;}
        ~safeBoolBase() {}
};


/*---------------------------------------------------------------------------*\
                        Class safeBool Declaration
\*---------------------------------------------------------------------------*/

template <typename T=void>
class safeBool
:
    private safeBoolBase
{
  // private or protected inheritance is very important here as it triggers the
  // access control violation in main.

public:

    operator boolType() const
    {
        return (static_cast<const T*>(this))->test()
            ? &safeBoolBase::thisTypeDoesNotSupportComparisons : 0;
    }

protected:

    //- Destructor
    ~safeBool()
    {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "safeBoolI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
