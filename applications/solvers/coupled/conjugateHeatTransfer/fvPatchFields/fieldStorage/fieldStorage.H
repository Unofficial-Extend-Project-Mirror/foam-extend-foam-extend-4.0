/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    fieldStorage

Description
    Region couple patch field

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

SourceFiles
    fieldStorage.C

\*---------------------------------------------------------------------------*/

#ifndef fieldStorage_H
#define fieldStorage_H

#include "coupledFvPatchField.H"
#include "fvPatchFieldMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class fieldStorage Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class fieldStorage
{
    // Private data

        // Reference to time
        const runTime& time_;

        //- original patch field
        Field<Type> storedField_;

        //- Current time index used to store originalPatchField_
        label curTimeIndex_;


public:

    // Constructors

        //- Construct from patch and internal field
        fieldStorage();

    // Member functions

        //- Store field
        void store(const Field<Type>&) const;

        //- Return stored field
        const Field<Type>& originalPatchField() const;

        //- Map (and resize as needed) from self given a mapping object
        void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        void rmap
        (
            const fvPatchField<Type>&,
            const labelList&
        );

        //- Clear storage
        void clear() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "fieldStorage.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
