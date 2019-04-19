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

#ifndef doubleFloat_H
#define doubleFloat_H

#include "label.H"
#include "products.H"

#include <cmath>

#if defined(darwin) && defined(__clang__)
#ifndef DUMMY_SCALAR_FUNCTIONS
#define DUMMY_SCALAR_FUNCTIONS
inline float j0f(float x) { return float(j0(double(x))); }
inline float j1f(float x) { return float(j1(double(x))); }
inline float y0f(float x) { return float(y0(double(x))); }
inline float y1f(float x) { return float(y1(double(x))); }
inline float jnf(const int n, const float s) { return float(jn(n, double(s))); }
inline float ynf(const int n, const float s) { return float(yn(n, double(s))); }

inline long double j0l(float x) { return double(j0(double(x))); }
inline long double j1l(float x) { return double(j1(double(x))); }
inline long double y0l(float x) { return double(y0(double(x))); }
inline long double y1l(float x) { return double(y1(double(x))); }
inline long double jnl(const int n, const float s) { return double(jn(n, double(s))); }
inline long double ynl(const int n, const float s) { return double(yn(n, double(s))); }
#endif
#endif // darwin

#if defined(mingw)
#ifndef DUMMY_SCALAR_FUNCTIONS
#define DUMMY_SCALAR_FUNCTIONS
inline float j0f(float x) { return _j0(x); }
inline float j1f(float x) { return _j1(x); }
inline float y0f(float x) { return _y0(x); }
inline float y1f(float x) { return _y1(x); }
inline float jnf(const int n, float s) { return _jn(n, s); }
inline float ynf(const int n, float s) { return _yn(n, s); }

inline double j0(double x) { return _j0(x); }
inline double j1(double x) { return _j1(x); }
inline double y0(double x) { return _y0(x); }
inline double y1(double x) { return _y1(x); }
inline double jn(const int n, double s) { return _jn(n, s); }
inline double yn(const int n, double s) { return _yn(n, s); }

inline long double j0l(long double x) { return _j0(x); }
inline long double j1l(long double x) { return _j1(x); }
inline long double y0l(long double x) { return _y0(x); }
inline long double y1l(long double x) { return _y1(x); }
inline long double jnl(const int n, long double s) { return _jn(n, s); }
inline long double ynl(const int n, long double s) { return _yn(n, s); }
#endif
#endif // mingw

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Cmpt>
class typeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


template<class Cmpt>
class symmTypeOfRank<Cmpt, 0>
{
public:

    typedef Cmpt type;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
inline bool equal(const T& s1, const T& s2)
{
    return s1 == s2;
}


#define MAXMINPOW(retType, type1, type2)          \
                                                  \
MAXMIN(retType, type1, type2)                     \
                                                  \
inline retType pow(const type1 s, const type2 e)  \
{                                                 \
    return ::pow(s, e);                           \
}


MAXMINPOW(long double, long double, long double)
MAXMINPOW(long double, long double, double)
MAXMINPOW(long double, long double, float)
MAXMINPOW(long double, double, long double)
MAXMINPOW(long double, float, long double)
MAXMINPOW(long double, long double, int)
MAXMINPOW(long double, int, long double)
MAXMINPOW(long double, long double, long)
MAXMINPOW(long double, long, long double)

MAXMINPOW(double, double, double)
MAXMINPOW(double, double, float)
MAXMINPOW(double, float, double)
MAXMINPOW(double, double, int)
MAXMINPOW(double, int, double)
MAXMINPOW(double, double, long)
MAXMINPOW(double, long, double)

MAXMINPOW(float, float, float)
MAXMINPOW(float, float, int)
MAXMINPOW(float, int, float)
MAXMINPOW(float, float, long)
MAXMINPOW(float, long, float)

#undef MAXMINPOW


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
