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
    Foam::Hash

Description
    Hash function class for primitives.  All non-primitives used to hash
    entries on hash tables likely need a specialized version of this class.

\*---------------------------------------------------------------------------*/

#ifndef Hash_H
#define Hash_H

#include <inttypes.h>

#include "label.H"
#include "uLabel.H"
#include "Hasher.H"
#include "pTraits.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                            Class Hash Declaration
\*---------------------------------------------------------------------------*/

template<class PrimitiveType>
class Hash
{
public:

    Hash()
    {}

    unsigned operator()(const PrimitiveType& p, unsigned seed) const
    {
        return Hasher(&p, sizeof(p), seed);
    }

    unsigned operator()(const PrimitiveType& p) const
    {
        return Hasher(&p, sizeof(p));
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Hash specialization for hashing labels
template<>
class Hash<label>
{
public:

    Hash()
    {}

    //- Incrementally hash a label.
    //  This will necessarily return a different value than the
    //  non-incremental version.
    unsigned operator()(const label& p, unsigned seed) const
    {
        return Hasher(&p, sizeof(label), seed);
    }

    //- Return the unsigned representation of a label.
    //  This helps if people have relied on the hash value corresponding to
    //  the natural order.
    unsigned operator()(const label& p) const
    {
        return p;
    }
};


//- Hash specialization for hashing pointer addresses.
//  Treat a pointer like a intptr_t.
//  This should work for both 32-bit and 64-bit pointers.
template<>
class Hash<void*>
{
public:

    Hash()
    {}

    unsigned operator()(const void* const& p, unsigned seed) const
    {
        return Hash<intptr_t>()(intptr_t(p), seed);
    }

    unsigned operator()(const void* const& p) const
    {
        return Hash<intptr_t>()(intptr_t(p));
    }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
