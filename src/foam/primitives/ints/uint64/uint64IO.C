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

#include "error.H"

#include "uint64.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const uint64_t val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, uint64_t& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        i = uint64_t(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, uint64_t&)", is)
            << "wrong token type - expected uint64_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, uint64_t&)");

    return is;
}


uint64_t Foam::readUint64(Istream& is)
{
    uint64_t val;
    is >> val;

    return val;
}


bool Foam::read(const char* buf, uint64_t& s)
{
    char *endptr = nullptr;
    long l = strtol(buf, &endptr, 10);
    s = uint64_t(l);
    return (*endptr == 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const uint64_t i)
{
    os.write(label(i));
    os.check("Ostream& operator<<(Ostream&, const uint64_t)");
    return os;
}

#if WM_ARCH_OPTION == 64 && darwin && __clang__
Foam::Ostream& Foam::operator<<(Ostream& os, const uintptr_t i)
{
    os.write(label(i));
    os.check("Ostream& operator<<(Ostream&, const uintptr_t)");
    return os;
}
#endif

// ************************************************************************* //
