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

#include "writeRegisteredObject.H"
#include "dictionary.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(writeRegisteredObject, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::writeRegisteredObject::writeRegisteredObject
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    const bool loadFromFiles
)
:
    name_(name),
    exclusiveWriting_(false),
    obr_(obr),
    objectNames_()
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::writeRegisteredObject::~writeRegisteredObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::writeRegisteredObject::read(const dictionary& dict)
{
    dict.lookup("objectNames") >> objectNames_;
    dict.readIfPresent("exclusiveWriting", exclusiveWriting_);
}


void Foam::writeRegisteredObject::execute()
{
    // Do nothing - only valid on write
}


void Foam::writeRegisteredObject::end()
{
    // Do nothing - only valid on write
}


void Foam::writeRegisteredObject::timeSet()
{
    // Do nothing - only valid on write
}


void Foam::writeRegisteredObject::write()
{
    Info<< type() << " " << name_ << " output:" << nl;

    DynamicList<word> allNames(obr_.toc().size());
    forAll(objectNames_, i)
    {
        wordList names(obr_.names<regIOobject>(objectNames_[i]));

        if (names.size())
        {
            allNames.append(names);
        }
        else
        {
            WarningIn("Foam::writeRegisteredObject::write()")
                << "Object " << objectNames_[i] << " not found in "
                << "database. Available objects:" << nl << obr_.sortedToc()
                << endl;
        }
    }

    forAll(allNames, i)
    {
        regIOobject& obj =
            const_cast<regIOobject&>
            (
                obr_.lookupObject<regIOobject>(allNames[i])
            );

        if (exclusiveWriting_)
        {
            // Switch off automatic writing to prevent double write
            obj.writeOpt() = IOobject::NO_WRITE;
        }

        Info<< "    writing object " << obj.name() << nl << endl;

        obj.write();
    }
}


// ************************************************************************* //
