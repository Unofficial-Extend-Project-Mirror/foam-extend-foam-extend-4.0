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
    Foam::Constant

Description
    Templated basic entry that holds a constant value.

    Usage - for entry <entryName> having the value <value>:
    @verbatim
        <entryName>    constant  <value>
    @endverbatim

SourceFiles
    Constant.C

\*---------------------------------------------------------------------------*/

#ifndef Constant_H
#define Constant_H

#include "DataEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class Constant;

template<class Type>
Ostream& operator<<
(
    Ostream&,
    const Constant<Type>&
);

/*---------------------------------------------------------------------------*\
                           Class Constant Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class Constant
:
    public DataEntry<Type>
{
    // Private data

        //- Constant value
        Type value_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const Constant<Type>&);


public:

    // Runtime type information
    TypeName("constant");


    // Constructors

        //- Construct from entry name and Istream
        Constant(const word& entryName, Istream& is);

        //- Copy constructor
        Constant(const Constant<Type>& cnst);


    //- Destructor
    virtual ~Constant();


    // Member Functions

        //- Return constant value
        Type value(const scalar) const;

        //- Integrate between two values
        Type integrate(const scalar x1, const scalar x2) const;


    //- Ostream Operator
    friend Ostream& operator<< <Type>
    (
        Ostream&,
        const Constant<Type>&
    );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
Constant<label>::Constant(const word& entryName, Istream& is);

template<>
Constant<scalar>::Constant(const word& entryName, Istream& is);


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "Constant.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
