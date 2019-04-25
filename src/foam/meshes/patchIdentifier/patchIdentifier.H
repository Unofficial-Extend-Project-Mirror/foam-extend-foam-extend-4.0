/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::patchIdentifier

Description
    Identifies patch by name, patch index and physical type

SourceFiles
    patchIdentifier.C

\*---------------------------------------------------------------------------*/

#ifndef patchIdentifier_H
#define patchIdentifier_H

#include "wordList.H"
#include "label.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;

// Forward declaration of friend functions and operators
class patchIdentifier;

Ostream& operator<<(Ostream&, const patchIdentifier&);


/*---------------------------------------------------------------------------*\
                           Class patchIdentifier Declaration
\*---------------------------------------------------------------------------*/

class patchIdentifier
{
    // Private data

        //- Name of patch
        word name_;

        //- Index of patch in boundary
        label index_;

        //- Optional physical type
        mutable word physicalType_;

        //- Optional groups patch belongs to
        wordList inGroups_;

public:

    // Constructors

        //- Construct from components
        patchIdentifier
        (
            const word& name,
            const label index,
            const word& physicalType = word::null,
            const wordList& inGroups = wordList()
        );

        //- Construct from dictionary
        patchIdentifier
        (
            const word& name,
            const dictionary&,
            const label index
        );

        //- Construct from geometric patch, resetting the index
        patchIdentifier
        (
            const patchIdentifier&,
            const label index
        );


    //- Destructor
    virtual ~patchIdentifier();


    // Member Functions

        //- Return name
        const word& name() const
        {
            return name_;
        }

        //- Return name for modification
        word& name()
        {
            return name_;
        }

        //- Return the optional physical type of the patch
        const word& physicalType() const
        {
            return physicalType_;
        }

        //- Return the optional physical type of the patch for modification
        word& physicalType()
        {
            return physicalType_;
        }

        //- Return the index of this patch in the boundaryMesh
        label index() const
        {
            return index_;
        }

        //- Return the index of this patch in the boundaryMesh for modification
        label& index()
        {
            return index_;
        }

        //- Return the optional groups patch belongs to
        const wordList& inGroups() const
        {
            return inGroups_;
        }

        //- Return the optional groups patch belongs to for modification
        wordList& inGroups()
        {
            return inGroups_;
        }

        //- Test if in group
        bool inGroup(const word&) const;

        //- Write patchIdentifier as a dictionary
        void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const patchIdentifier&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
