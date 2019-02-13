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

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "objectRefinement.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(objectRefinement, 0);
defineRunTimeSelectionTable(objectRefinement, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

objectRefinement::objectRefinement()
:
    name_(),
    cellSize_(-1.0),
    additionalRefLevel_(0),
    refThickness_(0.0)
{}

objectRefinement::objectRefinement
(
    const word& name,
    const dictionary& dict
)
:
    name_(name),
    cellSize_(-1.0),
    additionalRefLevel_(0),
    refThickness_(0.0)
{
    if( dict.found("cellSize") )
    {
        cellSize_ = readScalar(dict.lookup("cellSize"));

        if( cellSize_ < 0.0 )
        {
            FatalErrorIn
            (
                "objectRefinement::objectRefinement"
                "(const word&, const dictionary&)"
            ) << "Specified cell size for object " << name_
              << " is negative" << exit(FatalError);
        }
    }
    else if( dict.found("additionalRefinementLevels") )
    {
        additionalRefLevel_ =
            readLabel(dict.lookup("additionalRefinementLevels"));

        if( additionalRefLevel_ < 0 )
        {
            FatalErrorIn
            (
                "objectRefinement::objectRefinement"
                "(const word&, const dictionary&)"
            ) << "Specified additionalRefinementLevel for object " << name_
              << " is negative" << exit(FatalError);
        }
    }

    if( dict.found("refinementThickness") )
        refThickness_ = readScalar(dict.lookup("refinementThickness"));
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

objectRefinement::~objectRefinement()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const objectRefinement& obr)
{
    os << obr.name() << nl;
    obr.writeDict(os, true);
    return os;
}

void objectRefinement::calculateAdditionalRefLevels(const scalar globalCellSize)
{
    if( cellSize_ < 0.0 || additionalRefLevel_ != 0 )
        return;

    scalar s = globalCellSize;

    label nMarked;
    do
    {
        nMarked = 0;

        if( cellSize_ <= s * (1.+SMALL) )
        {
            ++nMarked;
            ++additionalRefLevel_;
        }

        s /= 2.0;

    } while( nMarked != 0 );

}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
