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

#include "fieldAverageItem.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

template<>
const char* Foam::NamedEnum<Foam::fieldAverageItem::baseType, 2>::names[] =
{
    "iteration",
    "time"
};

const Foam::NamedEnum<Foam::fieldAverageItem::baseType, 2>
    Foam::fieldAverageItem::baseTypeNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fieldAverageItem::fieldAverageItem()
:
    fieldName_("unknown"),
    mean_(0),
    prime2Mean_(0),
    base_(ITER)
{}


Foam::fieldAverageItem::fieldAverageItem(const fieldAverageItem& faItem)
:
    fieldName_(faItem.fieldName_),
    mean_(faItem.mean_),
    prime2Mean_(faItem.prime2Mean_),
    base_(faItem.base_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fieldAverageItem::~fieldAverageItem()
{}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::fieldAverageItem::operator=(const fieldAverageItem& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "Foam::fieldAverageItem::operator=(const Foam::fieldAverageItem&)"
        )   << "Attempted assignment to self" << nl
            << abort(FatalError);
    }

    // Set updated values
    fieldName_ = rhs.fieldName_;
    mean_ = rhs.mean_;
    prime2Mean_ = rhs.prime2Mean_;
    base_ = rhs.base_;
}


// ************************************************************************* //
