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
    Foam::morphFieldMapper

Description
    Abstract base class to hold the Field mapping for mesh morphs.

\*---------------------------------------------------------------------------*/

#ifndef morphFieldMapper_H
#define morphFieldMapper_H

#include "Field.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class morphFieldMapper Declaration
\*---------------------------------------------------------------------------*/

class morphFieldMapper
:
    public FieldMapper
{

public:

    // Constructors

        //- Null constructor
        morphFieldMapper()
        {}


    // Destructor

        virtual ~morphFieldMapper()
        {}


    // Member Functions

        //- Return size of field before mapping
        virtual label sizeBeforeMapping() const = 0;

        //- Are there any inserted objects
        virtual bool insertedObjects() const = 0;

        //- Return list of inserted objects
        virtual const labelList& insertedObjectLabels() const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
