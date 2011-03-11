/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-6 H. Jasak All rights reserved
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
    Generic field type.

\*---------------------------------------------------------------------------*/

#include "FadOne.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * */

// Return a string representation of a FadOne
template<int nVars>
word name(const FadOne<nVars>&);


// Allow binary writing for arrays of complex
template<int nVars>
inline bool writeBinary(const FadOne<nVars>*)
{
    return true;
}


template<int nVars>
inline const FadOne<nVars>& max
(
    const FadOne<nVars>&,
    const FadOne<nVars>&
);

template<int nVars>
inline const FadOne<nVars>& min
(
    const FadOne<nVars>&,
    const FadOne<nVars>&
);

template<int nVars>
inline FadOne<nVars> limit
(
    const FadOne<nVars>&,
    const FadOne<nVars>&
);


template<int nVars>
inline FadOne<nVars>& setComponent(FadOne<nVars>& u, const direction)
{
    return u;
}

template<int nVars>
inline FadOne<nVars> component(const FadOne<nVars>& u, const direction)
{
    return u;
}

template<int nVars>
inline FadOne<nVars> mag(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> sin(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> asin(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> cos(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> acos(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> tan(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> atan(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> tanh(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> exp(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> log(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> sqr(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> sqrt(const FadOne<nVars>& u);

template<int nVars>
inline FadOne<nVars> pow
(
    const FadOne<nVars>&,
    const FadOne<nVars>&
);

template<int nVars>
inline FadOne<nVars> pow
(
    const FadOne<nVars>&,
    const double
);

template<int nVars>
inline FadOne<nVars> pow
(
    const double&,
    const FadOne<nVars>&
);

template<int nVars>
inline label sign(const FadOne<nVars>& u);

template<int nVars>
inline label pos(const FadOne<nVars>& u);

template<int nVars>
inline label neg(const FadOne<nVars>& u);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
