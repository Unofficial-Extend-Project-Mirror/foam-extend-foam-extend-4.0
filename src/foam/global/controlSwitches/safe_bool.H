/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
    Foam::safe_bool

Description
    Implementation of the Safe Bool Idiom.
    This is the version without virtual functions.

    Author: Bjorn Karlsson
    http://www.artima.com/cppsource/safebool.html

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef safe_bool_H
#define safe_bool_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class safe_bool_base Declaration
\*---------------------------------------------------------------------------*/
class safe_bool_base
{
    public:
        typedef void (safe_bool_base::*bool_type)() const;
        void this_type_does_not_support_comparisons() const {}

    protected:

        safe_bool_base() {}
        safe_bool_base(const safe_bool_base&) {}
        safe_bool_base& operator=(const safe_bool_base&) {return *this;}
        ~safe_bool_base() {}
};

/*---------------------------------------------------------------------------*\
                          Class safe_bool Declaration
\*---------------------------------------------------------------------------*/

template <typename T=void>
class safe_bool
:
    private safe_bool_base
{
  // private or protected inheritance is very important here as it triggers the
  // access control violation in main.

public:
    operator bool_type() const
    {
	return (static_cast<const T*>(this))->boolean_test()
	    ? &safe_bool_base::this_type_does_not_support_comparisons : 0;
    }

protected:
    //- Destructor
    ~safe_bool()
    {}
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "safe_boolI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
