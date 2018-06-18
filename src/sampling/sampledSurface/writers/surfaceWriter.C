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

\*---------------------------------------------------------------------------*/

#include "surfaceWriter.H"

#include "MeshedSurfaceProxy.H"
#include "proxySurfaceWriter.H"

#include "HashTable.H"
#include "word.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(surfaceWriter, 0);
    defineRunTimeSelectionTable(surfaceWriter, word);
    defineRunTimeSelectionTable(surfaceWriter, wordDict);
    addNamedToRunTimeSelectionTable
    (
        surfaceWriter,
        surfaceWriter,
        word,
        null
    );
}


// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType)
{
    wordConstructorTable::iterator cstrIter =
        wordConstructorTablePtr_->find(writeType);

    if (cstrIter == wordConstructorTablePtr_->end())
    {
        if (MeshedSurfaceProxy<face>::canWriteType(writeType))
        {
            // generally unknown, but can be written via MeshedSurfaceProxy
            // use 'proxy' handler instead
            return autoPtr<surfaceWriter>(new proxySurfaceWriter(writeType));
        }

        if (cstrIter == wordConstructorTablePtr_->end())
        {
            FatalErrorIn
            (
                "surfaceWriter::New(const word&)"
            )   << "Unknown write type \"" << writeType << "\"\n\n"
                << "Valid write types : "
                << wordConstructorTablePtr_->sortedToc() << nl
                << "Valid proxy types : "
                << MeshedSurfaceProxy<face>::writeTypes() << endl
                << exit(FatalError);
        }
    }

    return autoPtr<surfaceWriter>(cstrIter()());
}


Foam::autoPtr<Foam::surfaceWriter>
Foam::surfaceWriter::New(const word& writeType, const dictionary& optDict)
{
    // find constructors with dictionary options
    wordDictConstructorTable::iterator cstrIter =
        wordDictConstructorTablePtr_->find(writeType);

    if (cstrIter == wordDictConstructorTablePtr_->end())
    {
        // revert to versions without options
        return Foam::surfaceWriter::New(writeType);
    }

    return autoPtr<surfaceWriter>(cstrIter()(optDict));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfaceWriter::surfaceWriter()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceWriter::~surfaceWriter()
{}


// ************************************************************************* //
