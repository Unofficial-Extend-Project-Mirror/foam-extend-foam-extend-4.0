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

#include "dlLibraryTable.H"
#include "OSspecific.H"
#include "int.H"
#include <dlfcn.h> 

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dlLibraryTable, 0);
}

Foam::dlLibraryTable Foam::dlLibraryTable::loadedLibraries;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
{}


Foam::dlLibraryTable::dlLibraryTable
(
    const dictionary& dict,
    const word& libsEntry
)
{
    open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
    forAllReverse(libPtrs_, i)
    {
        if (libPtrs_[i])
        {
            dlClose(libPtrs_[i]);
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open
(
    const fileName& functionLibName,
    const bool verbose
)
{
    if (functionLibName.size())
    {
        void* functionLibPtr =
            dlopen(functionLibName.c_str(), RTLD_LAZY|RTLD_GLOBAL);

#ifdef darwin
        // If failing to load under OS X, let's try some obvious variations
        // before giving up completely
        fileName osxFileName(functionLibName);

        if(!functionLibPtr && functionLibName.ext()=="so")
        {
            osxFileName=functionLibName.lessExt()+".dylib";

            functionLibPtr =
                dlopen(osxFileName.c_str(), RTLD_LAZY|RTLD_GLOBAL);
        }

        // If unsuccessful, which might be the case under Mac OSX 10.11 (El
        // Capitan) with System Integrity Protection (SIP) enabled, let's try
        // building a full path using well-known environment variables. This is
        // the last resort, unless you provide the full pathname yourself.
        if (!functionLibPtr)
        {
            fileName l_LIBBIN_Name =
                getEnv("FOAM_LIBBIN")/osxFileName;
            functionLibPtr =
                dlopen(l_LIBBIN_Name.c_str(), RTLD_LAZY|RTLD_GLOBAL);
        }
        if (!functionLibPtr)
        {
            fileName l_SITE_LIBBIN_Name =
                getEnv("FOAM_SITE_LIBBIN")/osxFileName;
            functionLibPtr =
                dlopen(l_SITE_LIBBIN_Name.c_str(), RTLD_LAZY|RTLD_GLOBAL);
        }
        if (!functionLibPtr)
        {
            fileName l_USER_LIBBIN_Name =
                getEnv("FOAM_USER_LIBBIN")/osxFileName;
            functionLibPtr =
                dlopen(l_USER_LIBBIN_Name.c_str(), RTLD_LAZY|RTLD_GLOBAL);
        }
#elif defined mingw
        if(!functionLibPtr && functionLibName.ext()=="so") {
            fileName lName=functionLibName.lessExt()+".dll";
            functionLibPtr =
                dlopen(lName.c_str(), RTLD_LAZY|RTLD_GLOBAL);
        }
#endif
        if (!functionLibPtr)
        {
            WarningIn
            (
                "dlLibraryTable::open(const fileName& functionLibName)"
            )   << "could not load " << dlerror()
                << endl;

            return false;
        }
        else
        {
            if (!loadedLibraries.found(functionLibPtr))
            {
                loadedLibraries.insert(functionLibPtr, functionLibName);
                return true;
            }
            else
            {
                return false;
            }
        }
    }
    else
    {
        return false;
    }
}


bool Foam::dlLibraryTable::close
(
    const fileName& functionLibName,
    const bool verbose
)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libNames_[i] == functionLibName)
        {
            index = i;
            break;
        }
    }

    if (index != -1)
    {
        bool ok = dlClose(libPtrs_[index]);

        libPtrs_[index] = NULL;
        libNames_[index] = fileName::null;

        if (!ok)
        {
            if (verbose)
            {
                WarningInFunction
                    << "could not close " << functionLibName
                    << endl;
            }

            return false;
        }

        return true;
    }
    return false;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& functionLibName)
{
    label index = -1;
    forAllReverse(libNames_, i)
    {
        if (libNames_[i] == functionLibName)
        {
            index = i;
            break;
        }
    }

    if (index != -1)
    {
        return libPtrs_[index];
    }
    return NULL;
}


bool Foam::dlLibraryTable::open
(
    const dictionary& dict,
    const word& libsEntry
)
{
    if (dict.found(libsEntry))
    {
        fileNameList libNames(dict.lookup(libsEntry));

        bool allOpened = !libNames.empty();

        forAll(libNames, i)
        {
            allOpened = dlLibraryTable::open(libNames[i]) && allOpened;
        }

        return allOpened;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
