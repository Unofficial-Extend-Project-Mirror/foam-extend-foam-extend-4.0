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

#include "prefixOSstream.H"
#include "Pstream.H"
#include "token.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

prefixOSstream::prefixOSstream
(
    ostream& os,
    const string& name,
    streamFormat format,
    versionNumber version,
    compressionType compression
)
:
    OSstream(os, name, format, version, compression),
    printPrefix_(true),
    prefix_("")
{}


inline void prefixOSstream::checkWritePrefix()
{
    if (printPrefix_ && prefix_.size())
    {
        OSstream::write(prefix_.c_str());
        printPrefix_ = false;
    }
}


Ostream& prefixOSstream::write(const token&)
{
    return *this;
}


Ostream& prefixOSstream::write(const char c)
{
    checkWritePrefix();
    OSstream::write(c);

    if (c == token::NL)
    {
        printPrefix_ = true;
    }

    return *this;
}


Ostream& prefixOSstream::write(const char* s)
{
    checkWritePrefix();
    OSstream::write(s);

    size_t sl = strlen(s);
    if (sl && s[sl - 1] == token::NL)
    {
        printPrefix_ = true;
    }

    return *this;
}


Ostream& prefixOSstream::write(const word& w)
{
    checkWritePrefix();
    return OSstream::write(w);
}


Ostream& prefixOSstream::write(const string& s)
{
    checkWritePrefix();
    return OSstream::write(s);
}


Ostream& prefixOSstream::write(const label l)
{
    checkWritePrefix();
    return OSstream::write(l);
}


Ostream& prefixOSstream::write(const floatScalar s)
{
    checkWritePrefix();
    return OSstream::write(s);
}


Ostream& prefixOSstream::write(const doubleScalar s)
{
    checkWritePrefix();
    return OSstream::write(s);
}


Ostream& prefixOSstream::write(const char* buf, std::streamsize count)
{
    checkWritePrefix();
    return OSstream::write(buf, count);
}


//- Add indentation characters
void prefixOSstream::indent()
{
    checkWritePrefix();
    OSstream::indent();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
