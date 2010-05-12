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

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

IOstream::streamFormat IOstream::formatEnum(const word& format)
{
    if (format == "ascii")
    {
        return IOstream::ASCII;
    }
    else if (format == "binary")
    {
        return IOstream::BINARY;
    }
    else
    {
        WarningIn("IOstream::formatEnum(const word&)")
            << "bad format specifier "
            << format << " using ASCII"
            << endl;

        return IOstream::ASCII;
    }
}


IOstream::compressionType IOstream::compressionEnum(const word& compression)
{
    if (compression == "uncompressed")
    {
        return IOstream::UNCOMPRESSED;
    }
    else if (compression == "compressed")
    {
        return IOstream::COMPRESSED;
    }
    else
    {
        WarningIn("IOstream::compressionEnum(const word&)")
            << "bad compression specifier "
            << '\'' << compression << '\''
            << ", use 'compressed' or 'uncompressed'.  "
               "Defaulting to uncompressed"
            << endl;

        return IOstream::UNCOMPRESSED;
    }
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const IOstream::streamFormat& sf)
{
    if (sf == IOstream::ASCII)
    {
        os << "ascii";
    }
    else
    {
        os << "binary";
    }

    return os;
}


#if defined (__GNUC__)
template<>
#endif
Ostream& operator<<(Ostream& os, const InfoProxy<IOstream>& iostr)
{
    iostr.t_.print(os);
    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
