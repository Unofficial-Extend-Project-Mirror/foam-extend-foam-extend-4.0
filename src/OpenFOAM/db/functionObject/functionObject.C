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

\*---------------------------------------------------------------------------*/

#include "functionObject.H"
#include "dictionary.H"
#include "dlLibraryTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineRunTimeSelectionTable(functionObject, dictionary);
    int functionObject::debug(::Foam::debug::debugSwitch("functionObject", 0));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObject::functionObject()
{}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObject> Foam::functionObject::New
(
    const word& name,
    const Time& t,
    const dictionary& functionDict
)
{
    word functionType(functionDict.lookup("type"));

    if (debug)
    {
        Info<< "Selecting function " << functionType << endl;
    }

    dlLibraryTable::open
    (
        functionDict,
        "functionObjectLibs",
        dictionaryConstructorTablePtr_
    );

    if (!dictionaryConstructorTablePtr_)
    {
        FatalErrorIn
        (
            "functionObject::New"
            "(const word& functionType, const Time&, const dictionary&)"
        )   << "Unknown function type "
            << functionType << endl << endl
            << "Table of functionObjects is empty"
            << exit(FatalError);
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(functionType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "functionObject::New"
            "(const word& functionType, const Time&, const dictionary&)"
        )   << "Unknown function type "
            << functionType << endl << endl
            << "Valid functions are : " << endl
            << dictionaryConstructorTablePtr_->toc()
            << exit(FatalError);
    }

    return autoPtr<functionObject>(cstrIter()(name, t, functionDict));
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObject::~functionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::functionObject> Foam::functionObject::iNew::operator()
(
    const word& name,
    Istream& is
) const
{
    dictionary dict(is);
    return functionObject::New(name, time_, dict);
}


// ************************************************************************* //
