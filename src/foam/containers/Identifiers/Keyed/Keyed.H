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
    Foam::Keyed

Description
    A container with an integer key attached to any item.

    The key can useful for sorting.

SourceFiles
    KeyedI.H

\*---------------------------------------------------------------------------*/

#ifndef Keyed_H
#define Keyed_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class T> class Keyed;

template<class T> Istream& operator>>(Istream&, Keyed<T>&);
template<class T> Ostream& operator<<(Ostream&, const Keyed<T>&);

/*---------------------------------------------------------------------------*\
                       Class Keyed Declaration
\*---------------------------------------------------------------------------*/

template<class T>
class Keyed
:
    public T
{
    // Private data

        label key_;

public:

    // Static Members

        //- Add labels to a list of values
        inline static List<Keyed<T> > createList
        (
            const List<T>&,
            const label key=0
        );

        //- Add labels to a list of values
        inline static List<Keyed<T> > createList
        (
            const List<T>&,
            const List<label>& keys
        );


    // Constructors

        //- Construct null
        inline Keyed();

        //- Construct as a copy of item, with a key
        inline Keyed(const T& item, const label key=0);

        //- Construct by transferring the item, with a key
        inline Keyed(const Xfer<T>& item, const label key=0);

        //- Construct from Istream
        inline Keyed(Istream&);


    // Member Functions

        // Access

            //- Return const access to the integer key
            inline label key() const;

            //- Return non-const access to the integer key
            inline label& key();


    // IOstream Operators

        friend Istream& operator>> <T>(Istream&, Keyed<T>&);
        friend Ostream& operator<< <T>(Ostream&, const Keyed<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "KeyedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
