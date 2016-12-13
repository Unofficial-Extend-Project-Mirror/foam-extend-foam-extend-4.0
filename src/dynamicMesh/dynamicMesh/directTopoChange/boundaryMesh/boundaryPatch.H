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
    Foam::boundaryPatch

Description
    Like polyPatch but without reference to mesh. patchIdentifier::index
    is not used. Used in boundaryMesh to hold data on patches.

SourceFiles
    boundaryPatch.C

\*---------------------------------------------------------------------------*/

#ifndef boundaryPatch_H
#define boundaryPatch_H

#include "patchIdentifier.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class boundaryPatch Declaration
\*---------------------------------------------------------------------------*/

class boundaryPatch
:
    public patchIdentifier
{
    // Private data

        label size_;
        label start_;

public:

    // Constructors

        //- Construct from components
        boundaryPatch
        (
            const word& name,
            const label index,
            const label size,
            const label start,
            const word& physicalType = word::null
        );

        //- Construct from dictionary
        boundaryPatch
        (
            const word& name,
            const dictionary& dict,
            const label index
        );

        //- Construct as copy
        boundaryPatch(const boundaryPatch&);

        //- Clone
        autoPtr<boundaryPatch> clone() const;


    // Destructor

        ~boundaryPatch();


    // Member Functions

        label size() const
        {
            return size_;
        }

        label& size()
        {
            return size_;
        }

        label start() const
        {
            return start_;
        }

        label& start()
        {
            return start_;
        }


        //- Write dictionary
        virtual void write(Ostream&) const;


    // Ostream Operator

        friend Ostream& operator<<(Ostream&, const boundaryPatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
