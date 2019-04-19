/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    boundaryPatch

Description
    Writes information for a given boundary patch

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryPatch_H
#define boundaryPatch_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "boundaryPatchBase.H"
#include "typeInfo.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class boundaryPatch
:   public boundaryPatchBase
{

public:

    //- Runtime type information
    TypeName("patch");

    // Constructors

        boundaryPatch
        (
            const word&,
            const word&,
            const label,
            const label
        );

        boundaryPatch(const word& name, const dictionary&);

    // Return clone of the object
        virtual autoPtr<boundaryPatchBase> clone() const
        {
            return autoPtr<boundaryPatchBase>
            (
                new boundaryPatch
                (
                    patchName(),
                    patchType(),
                    patchSize(),
                    patchStart()
                )
            );
        }

        virtual autoPtr<boundaryPatchBase> clone(const boundaryPatch& wp) const
        {
            return autoPtr<boundaryPatchBase>
            (
                new boundaryPatch
                (
                    wp.patchName(),
                    wp.patchType(),
                    wp.patchSize(),
                    wp.patchStart()
                )
            );
        }

        //- Return as dictionary of entries
        dictionary dict() const;

    // Write
        //- write to Ostream
        void write(Ostream&) const;

        //- Write dictionary
        void writeDict(Ostream&) const;

    // Member operators

        Ostream& operator<<(Ostream&) const;

        Istream& operator>>(Istream&);

        void operator=(const boundaryPatch&);

        bool operator!=(const boundaryPatch&) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
