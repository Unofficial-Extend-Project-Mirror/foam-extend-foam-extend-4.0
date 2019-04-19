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

#include "IOobject.H"
#include "foamTime.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::IOobject, 0);

// * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * * //

bool Foam::IOobject::fileNameComponents
(
    const fileName& path,
    fileName& instance,
    fileName& local,
    word& name
)
{
    instance.clear();
    local.clear();
    name.clear();

    // called with directory
    if (isDir(path))
    {
        WarningIn
        (
            "IOobject::fileNameComponents"
            "("
                "const fileName&, "
                "fileName&, "
                "fileName&, "
                "word&"
            ")"
        )   << " called with directory: " << path << endl;

        return false;
    }

    if (path.isAbsolute())
    {
        string::size_type last = path.rfind('/');
        instance = path.substr(0, last);

        // Check afterwards
        name.string::operator=(path.substr(last+1));
    }
    else
    {
        string::size_type first = path.find('/');

        if (first == string::npos)
        {
            // no '/' found - no instance or local

            // check afterwards
            name.string::operator=(path);
        }
        else
        {
            instance = path.substr(0, first);

            string::size_type last = path.rfind('/');
            if (last > first)
            {
                // with local
                local = path.substr(first+1, last-first-1);
            }

            // check afterwards
            name.string::operator=(path.substr(last+1));
        }
    }


    // Check for valid (and stripped) name, regardless of the debug level
    if (name.empty() || string::stripInvalid<word>(name))
    {
        WarningIn
        (
            "IOobject::fileNameComponents"
            "("
                "const fileName&, "
                "fileName&, "
                "fileName&, "
                "word&"
            ")"
        )
            << "has invalid word for name: \"" << name
            << "\"\nwhile processing path: " << path << endl;

        return false;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::IOobject::IOobject
(
    const word& name,
    const Time& time,
    readOption ro,
    writeOption wo
)
:
    name_(name),
    headerClassName_(typeName),
    note_(),
    instance_(""),
    local_(),
    db_(time),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(false),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(registry.mangleFileName(name)),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const word& name,
    const fileName& instance,
    const fileName& local,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(registry.mangleFileName(name)),
    headerClassName_(typeName),
    note_(),
    instance_(instance),
    local_(local),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


Foam::IOobject::IOobject
(
    const fileName& path,
    const objectRegistry& registry,
    readOption ro,
    writeOption wo,
    bool registerObject
)
:
    name_(),
    headerClassName_(typeName),
    note_(),
    instance_(),
    local_(),
    db_(registry),
    rOpt_(ro),
    wOpt_(wo),
    registerObject_(registerObject),
    objState_(GOOD)
{
    if (!fileNameComponents(path, instance_, local_, name_))
    {
        FatalErrorIn
        (
            "IOobject::IOobject"
            "("
                "const fileName&, "
                "const objectRegistry&, "
                "readOption, "
                "writeOption, "
                "bool"
            ")"
        )
            << " invalid path specification"
            << exit(FatalError);
    }

    if (objectRegistry::debug)
    {
        Info<< "Constructing IOobject called " << name_
            << " of type " << headerClassName_
            << endl;
    }
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

Foam::IOobject::~IOobject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::objectRegistry& Foam::IOobject::db() const
{
    return db_;
}


const Foam::Time& Foam::IOobject::time() const
{
    return db_.time();
}


const Foam::fileName& Foam::IOobject::caseName() const
{
    return time().caseName();
}


Foam::word Foam::IOobject::group() const
{
    word::size_type i = name_.find_last_of('.');

    if (i == word::npos || i == 0)
    {
        return word::null;
    }
    else
    {
        return name_.substr(i+1, word::npos);
    }
}


Foam::word Foam::IOobject::member() const
{
    word::size_type i = name_.find_last_of('.');

    if (i == word::npos || i == 0)
    {
        return name_;
    }
    else
    {
        return name_.substr(0, i);
    }
}


const Foam::fileName& Foam::IOobject::rootPath() const
{
    return time().rootPath();
}


Foam::fileName Foam::IOobject::path() const
{
    if (instance().isAbsolute())
    {
        return instance();
    }
    else
    {
        return rootPath()/caseName()/instance()/db_.dbDir()/local();
    }
}


Foam::fileName Foam::IOobject::path
(
    const word& instance,
    const fileName& local
) const
{
    // Note: can only be called with relative instance since is word type
    return rootPath()/caseName()/instance/db_.dbDir()/local;
}


Foam::fileName Foam::IOobject::filePath() const
{
    if (instance().isAbsolute())
    {
        fileName objectPath = instance()/name();
        if (isFile(objectPath))
        {
            return objectPath;
        }
        else
        {
            return fileName::null;
        }
    }
    else
    {
        fileName path = this->path();
        fileName objectPath = path/name();

        if (isFile(objectPath))
        {
            return objectPath;
        }
        else
        {
            if
            (
                time().processorCase()
             && (
                    instance() == time().system()
                 || instance() == time().constant()
                )
            )
            {
                fileName parentObjectPath =
                    rootPath()/time().globalCaseName()
                   /instance()/db_.dbDir()/local()/name();

                if (isFile(parentObjectPath))
                {
                    return parentObjectPath;
                }
            }

            if (!isDir(path))
            {
                word newInstancePath = time().findInstancePath
                (
                    instant(instance())
                );

                if (newInstancePath.size())
                {
                    fileName fName
                    (
                        rootPath()/caseName()
                       /newInstancePath/db_.dbDir()/local()/name()
                    );

                    if (isFile(fName))
                    {
                        return fName;
                    }
                }
            }
        }

        return fileName::null;
    }
}


Foam::Istream* Foam::IOobject::objectStream()
{
    return objectStream(filePath());
}


Foam::Istream* Foam::IOobject::objectStream(const fileName& fName)
{
    if (fName.size())
    {
        IFstream* isPtr = new IFstream(fName);

        if (isPtr->good())
        {
            return isPtr;
        }
        else
        {
            delete isPtr;
            return nullptr;
        }
    }
    else
    {
        return nullptr;
    }
}


bool Foam::IOobject::headerOk()
{
    bool ok = true;

    Istream* isPtr = objectStream();

    // If the stream has failed return
    if (!isPtr)
    {
        if (objectRegistry::debug)
        {
            Info
                << "IOobject::headerOk() : "
                << "file " << objectPath() << " could not be opened"
                << endl;
        }

        ok = false;
    }
    else
    {
        // Try reading header
        if (!readHeader(*isPtr))
        {
            if (objectRegistry::debug)
            {
                IOWarningIn("IOobject::headerOk()", (*isPtr))
                    << "failed to read header of file " << objectPath()
                    << endl;
            }

            ok = false;
        }
    }

    delete isPtr;

    return ok;
}


void Foam::IOobject::setBad(const string& s)
{
    if (objState_ != GOOD)
    {
        FatalErrorIn("IOobject::setBad(const string&)")
            << "recurrent failure for object " << s
            << exit(FatalError);
    }

    if (error::level)
    {
        Info<< "IOobject::setBad(const string&) : "
            << "broken object " << s << info() << endl;
    }

    objState_ = BAD;
}


void Foam::IOobject::operator=(const IOobject& io)
{
    name_ = io.name_;
    headerClassName_ = io.headerClassName_;
    note_ = io.note_;
    instance_ = io.instance_;
    local_ = io.local_;
    rOpt_ = io.rOpt_;
    wOpt_ = io.wOpt_;
    objState_ = io.objState_;
}


// ************************************************************************* //
