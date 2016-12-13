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
    Foam::polyPatchID

Description
    A class holds the data needed to identify a patch in a dynamic mesh.

    The patch is identified by name and its index in the boundary mesh
    is updated if the mesh has changed.

\*---------------------------------------------------------------------------*/

#ifndef polyPatchID_H
#define polyPatchID_H

#include "polyBoundaryMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class polyPatchID;
Ostream& operator<<(Ostream& os, const polyPatchID& p);


/*---------------------------------------------------------------------------*\
                           Class polyPatchID Declaration
\*---------------------------------------------------------------------------*/

class polyPatchID
{
    // Private data

        //- Patch name
        word name_;

        //- Patch index
        label index_;


public:

    // Constructors

        //- Construct from name
        polyPatchID
        (
            const word& name,
            const polyBoundaryMesh& bm,
            bool checkValid = false
        )
        :
            name_(name),
            index_(bm.findPatchID(name))
        {
            if (checkValid && index_ < 0)
            {
                FatalErrorIn
                (
                    "polyPatchID\n"
                    "(\n"
                    "    const word& name,\n"
                    "    const polyBoundaryMesh& bm,\n"
                    "    bool checkValid = false\n"
                    ")"
                )   << "Patch " << name_ << " not found.  Valid patch names: "
                    << bm.names()
                    << abort(FatalError);
            }
        }

        //- Construct from Istream
        polyPatchID
        (
            Istream& is,
            const polyBoundaryMesh& bm,
            bool checkValid = false
        )
        :
            name_(is),
            index_(bm.findPatchID(name_))
        {
            if (checkValid && index_ < 0)
            {
                FatalErrorIn
                (
                    "polyPatchID\n"
                    "(\n"
                    "    Istream& is\n"
                    "    const polyBoundaryMesh& bm,\n"
                    "    bool checkValid = false\n"
                    ")"
                )   << "Patch " << name_ << " not found.  Valid patch names: "
                    << bm.names()
                    << abort(FatalError);
            }
        }


    // Member Functions

        // Access

            //- Return name
            const word& name() const
            {
                return name_;
            }

            //- Return index
            label index() const
            {
                return index_;
            }

            //- Has the patch been found
            bool active() const
            {
                return index_ > -1;
            }


        // Edit

            //- Update
            void update(const polyBoundaryMesh& bm)
            {
                index_ = bm.findPatchID(name_);
            }


    // Ostream Operator

        friend Ostream& operator<<(Ostream& os, const polyPatchID& p)
        {
            os  << token::BEGIN_LIST
                << p.name_ << token::SPACE
                << p.index_
                << token::END_LIST;

            // Check state of Ostream
            os.check("Ostream& operator<<(Ostream&, const polyPatchID&)");

            return os;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
