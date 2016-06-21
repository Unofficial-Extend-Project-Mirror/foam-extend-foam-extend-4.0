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

\*---------------------------------------------------------------------------*/

#include "splitCell.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from cell number and parent
Foam::splitCell::splitCell(const label cellI, splitCell* parent)
:
    cellI_(cellI),
    parent_(parent),
    master_(NULL),
    slave_(NULL)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::splitCell::~splitCell()
{
    splitCell* myParent = parent();

    if (myParent)
    {
        // Make sure parent does not refer to me anymore.
        if (myParent->master() == this)
        {
            myParent->master() = NULL;
        }
        else if (myParent->slave() == this)
        {
            myParent->slave() = NULL;
        }
        else
        {
            FatalErrorIn("splitCell::~splitCell()") << "this not equal to"
                << " parent's master or slave pointer" << endl
                << "Cell:" << cellLabel() << abort(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::splitCell::isMaster() const
{
    splitCell* myParent = parent();

    if (!myParent)
    {
        FatalErrorIn("splitCell::isMaster()") << "parent not set"
            << "Cell:" << cellLabel() << abort(FatalError);

        return false;
    }
    else if (myParent->master() == this)
    {
        return true;
    }
    else if (myParent->slave() == this)
    {
        return false;
    }
    else
    {
        FatalErrorIn("splitCell::isMaster()") << "this not equal to"
            << " parent's master or slave pointer" << endl
            << "Cell:" << cellLabel() << abort(FatalError);

        return false;
    }
}


bool Foam::splitCell::isUnrefined() const
{
    return !master() && !slave();
}


Foam::splitCell* Foam::splitCell::getOther() const
{
    splitCell* myParent = parent();

    if (!myParent)
    {
        FatalErrorIn("splitCell::getOther()") << "parent not set"
            << "Cell:" << cellLabel() << abort(FatalError);

        return NULL;
    }
    else if (myParent->master() == this)
    {
        return myParent->slave();
    }
    else if (myParent->slave() == this)
    {
        return myParent->master();
    }
    else
    {
        FatalErrorIn("splitCell::getOther()") << "this not equal to"
            << " parent's master or slave pointer" << endl
            << "Cell:" << cellLabel() << abort(FatalError);

        return NULL;
    }
}


// ************************************************************************* //
