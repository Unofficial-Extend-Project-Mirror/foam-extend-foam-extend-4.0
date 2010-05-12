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

#include "functionObjectList.H"
#include "Time.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const bool execution
)
:
    HashPtrTable<functionObject>(),
    time_(t),
    foDict_(t.controlDict()),
    execution_(execution)
{}


Foam::functionObjectList::functionObjectList
(
    const Time& t,
    const dictionary& foDict,
    const bool execution
)
:
    HashPtrTable<functionObject>(),
    time_(t),
    foDict_(foDict),
    execution_(execution)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectList::~functionObjectList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjectList::start()
{
    if (execution_)
    {
        bool ok = false;

        if (foDict_.found("functions"))
        {
            HashPtrTable<functionObject> functions
            (
                foDict_.lookup("functions"),
                functionObject::iNew(time_)
            );

            transfer(functions);

            forAllIter(HashPtrTable<functionObject>, *this, iter)
            {
                ok = iter()->start() && ok;
            }
        }

        return ok;
    }
    else
    {
        return true;
    }
}


bool Foam::functionObjectList::execute()
{
    if (execution_)
    {
        bool ok = false;

        forAllIter(HashPtrTable<functionObject>, *this, iter)
        {
            ok = iter()->execute() && ok;
        }

        return ok;
    }
    else
    {
        return true;
    }
}


void Foam::functionObjectList::on()
{
    execution_ = true;
}


void Foam::functionObjectList::off()
{
    execution_ = false;
}


bool Foam::functionObjectList::read()
{
    bool read = false;

    if (foDict_.found("functions"))
    {
        HashPtrTable<dictionary> functionDicts(foDict_.lookup("functions"));

        // Update existing and add new functionObjects
        forAllConstIter(HashPtrTable<dictionary>, functionDicts, iter)
        {
            if (found(iter.key()))
            {
                read = find(iter.key())()->read(*iter()) && read;
            }
            else
            {
                functionObject* functionObjectPtr =
                    functionObject::New(iter.key(), time_, *iter()).ptr();

                functionObjectPtr->start();

                insert(iter.key(), functionObjectPtr);
            }
        }

        // Remove deleted functionObjects
        forAllIter(HashPtrTable<functionObject>, *this, iter)
        {
            if (!functionDicts.found(iter.key()))
            {
                erase(iter);
            }
        }
    }
    else
    {
        clear();
        read = true;
    }

    return read;
}


// ************************************************************************* //
