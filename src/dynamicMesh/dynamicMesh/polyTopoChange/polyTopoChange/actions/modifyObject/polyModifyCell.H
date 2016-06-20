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
    Foam::polyModifyCell

Description
    Class describing modification of a cell.

\*---------------------------------------------------------------------------*/

#ifndef polyModifyCell_H
#define polyModifyCell_H

#include "label.H"
#include "topoAction.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class polyModifyCell Declaration
\*---------------------------------------------------------------------------*/

class polyModifyCell
:
    public topoAction
{
    // Private data

        //- Cell ID
        label cellID_;

        //- Remove from current zone
        bool removeFromZone_;

        //- New zone ID
        label zoneID_;


public:

    // Static data members

        //- Runtime type information
        TypeName("modifyCell");


    // Constructors

        //- Construct null.  Used only for list construction
        polyModifyCell()
        :
            cellID_(-1),
            removeFromZone_(false),
            zoneID_(-1)
        {}

        //- Construct from components
        polyModifyCell
        (
            const label cellID,
            const bool removeFromZone,
            const label newZoneID
        )
        :
            cellID_(cellID),
            removeFromZone_(removeFromZone),
            zoneID_(newZoneID)
        {}

        //- Construct and return a clone
        virtual autoPtr<topoAction> clone() const
        {
            return autoPtr<topoAction>(new polyModifyCell(*this));
        }


    // Default Destructor

    // Member Functions

        //- Cell ID
        label cellID() const
        {
            return cellID_;
        }

        //- Does the cell belong to a zone?
        bool isInZone() const
        {
            return zoneID_ >= 0;
        }

        bool removeFromZone() const
        {
            return removeFromZone_;
        }

        //- Cell zone ID
        label zoneID() const
        {
            return zoneID_;
        }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
