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

Description
    Return the location of "directory" containing the file "name".
    Used in reading mesh data.

\*---------------------------------------------------------------------------*/

#include "Time.H"
#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::Time::findInstance
(
    const fileName& dir,
    const word& name,
    const IOobject::readOption rOpt
) const
{
    // Is the mesh data in the current time directory ?
    if
    (
        file(path()/timeName()/dir/name)
     && IOobject(name, timeName(), dir, *this).headerOk()
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance(const word& dir, const word& name) : "
                << "reading " << name
                << " from " << timeName()/dir
                << endl;
        }

        return timeName();
    }

    // Search back through the time directories list to find the time
    // closest to and lower than current time

    instantList ts = times();
    label i;

    for (i=ts.size()-1; i>=0; i--)
    {
        if (ts[i].value() <= timeOutputValue())
        {
            break;
        }
    }

    // Noting that the current directory has already been searched
    // for mesh data, start searching from the previously stored time directory

    if (i>=0)
    {
        for (label j=i; j>=0; j--)
        {
            if
            (
                file(path()/ts[j].name()/dir/name)
             && IOobject(name, ts[j].name(), dir, *this).headerOk()
            )
            {
                if (debug)
                {
                    Info<< "Time::findInstance(const word& dir, "
                        << "const word& name) : "
                        << "reading " << name
                        << " from " << ts[j].name()/dir
                        << endl;
                }

                return ts[j].name();
            }
        }
    }


    // If the mesh data is not in any of the time directories
    // Try in constant

    // Note. This needs to be a hard-coded constant, rather than the
    // constant function of the time, because the latter points to
    // the case constant directory in parallel cases

    if
    (
        file(path()/constant()/dir/name)
     && IOobject(name, constant(), dir, *this).headerOk()
    )
    {
        if (debug)
        {
            Info<< "Time::findInstance(const word& dir, "
                << "const word& name) : "
                << "reading " << name
                << " from " << constant()/dir
                << endl;
        }

        return constant();
    }

    if (rOpt == IOobject::MUST_READ)
    {
        FatalErrorIn("Time::findInstance(const word& dir, const word& name)")
            << "Cannot find file \"" << name << "\" in directory "
            << constant()/dir
            << exit(FatalError);
    }

    return constant();
}

// ************************************************************************* //
