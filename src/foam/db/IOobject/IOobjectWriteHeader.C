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

Description
    Writes the header description of the File to the stream
    associated with the File.

\*---------------------------------------------------------------------------*/

#include "IOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::IOobject::writeHeader(Ostream& os) const
{
    if (!os.good())
    {
        Info<< "IOobject::writeHeader(Ostream&) : "
            << "no stream open for write" << nl
            << os.info() << endl;

        return false;
    }

    writeBanner(os, false, true)
        << "FoamFile\n{\n"
        << "    version     " << os.version() << ";\n"
        << "    format      " << os.format() << ";\n"
        << "    class       " << type() << ";\n";

    if (note().size())
    {
        os  << "    note        " << note() << ";\n";
    }

    os  << "    location    " << instance()/local() << ";\n"
        << "    object      " << name() << ";\n"
        << "}" << nl;

    writeDivider(os) << endl;

    return true;
}


// ************************************************************************* //
