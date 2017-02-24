/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    sixDOFODE

Description
    Abstract base class for six-degrees-of-freedom (6DOF) ordinary differential
    equations

Author
    Dubravko Matijasevic, FSB Zagreb.  All rights reserved.
    Hrvoje Jasak, FSB Zagreb. All rights reserved.
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "objectRegistry.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::sixDOFODE> Foam::sixDOFODE::New(const IOobject& io)
{
    word sixDOFODETypeName;

    // Get object registry
    const objectRegistry& database = io.db();

    // Check whether the dictionary is in the database
    if (database.foundObject<IOdictionary>(io.name()))
    {
        sixDOFODETypeName =
            word
            (
                database.lookupObject<IOdictionary>(io.name()).lookup("type")
            );
    }
    else
    {
        sixDOFODETypeName = word(IOdictionary(io).lookup("type"));
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(sixDOFODETypeName);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "sixDOFODE::New(const IOobject& io)"
        )   << "Unknown sixDOFODE " << sixDOFODETypeName
            << endl << endl
            << "Valid sixDOFODE types are:" << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<sixDOFODE>(cstrIter()(io));
}


// ************************************************************************* //
