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

#include "int64.H"
#include "IOstreams.H"

#include <inttypes.h>
#include <sstream>
#include <cerrno>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const int64_t val)
{
    std::ostringstream buf;
    buf << val;
    return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, int64_t& i)
{
    token t(is);

    if (!t.good())
    {
        is.setBad();
        return is;
    }

    if (t.isLabel())
    {
        i = int64_t(t.labelToken());
    }
    else
    {
        is.setBad();
        FatalIOErrorIn("operator>>(Istream&, int64_t&)", is)
            << "wrong token type - expected int64_t, found " << t.info()
            << exit(FatalIOError);

        return is;
    }

    // Check state of Istream
    is.check("Istream& operator>>(Istream&, int64_t&)");

    return is;
}


int64_t Foam::readInt64(Istream& is)
{
    int64_t val;
    is >> val;

    return val;
}


bool Foam::read(const char* buf, int64_t& s)
{
    char *endptr = nullptr;
    errno = 0;
    intmax_t l = strtoimax(buf, &endptr, 10);
    s = int64_t(l);
    return (*endptr == 0) && (errno == 0);
}


Foam::Ostream& Foam::operator<<(Ostream& os, const int64_t i)
{
    os.write(label(i));
    os.check("Ostream& operator<<(Ostream&, const int64_t)");
    return os;
}

#if WM_ARCH_OPTION == 64 && darwin && __clang__
Foam::Istream& Foam::operator>>(Istream& is, long& i)
{
    return operator>>(is, reinterpret_cast<int64_t&>(i));
}

Foam::Ostream& Foam::operator<<(Ostream& os, const long i)
{
    os << int64_t(i);
    return os;
}
#endif

#if defined(mingw)
Foam::Istream& Foam::operator>>(Istream& is, off_t& i)
{
    return operator>>(is, i);
}

Foam::Ostream& Foam::operator<<(Ostream& os, const off_t i)
{
    os << i;
    return os;
}
#endif


// ************************************************************************* //
