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

Typedef
    Foam::scalar

Description
    Single floating point number identical to float or double depending on
    whether WM_SP or WM_DP is defined.

SourceFiles
    scalar.C

\*---------------------------------------------------------------------------*/

#ifndef scalar_H
#define scalar_H

#include "floatScalar.H"
#include "doubleScalar.H"
#include "longDoubleScalar.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#if defined(WM_SP)

// Define scalar as a float

namespace Foam
{
    typedef floatScalar scalar;

    static const scalar GREAT = floatScalarGREAT;
    static const scalar VGREAT = floatScalarVGREAT;
    static const scalar ROOTVGREAT = floatScalarROOTVGREAT;
    static const scalar SMALL = floatScalarSMALL;
    static const scalar VSMALL = floatScalarVSMALL;
    static const scalar ROOTVSMALL = floatScalarROOTVSMALL;

    scalar readScalar(Istream& is);
}

#elif defined(WM_DP)

// Define scalar as a double

namespace Foam
{
    typedef doubleScalar scalar;

    static const scalar GREAT = doubleScalarGREAT;
    static const scalar VGREAT = doubleScalarVGREAT;
    static const scalar ROOTVGREAT = doubleScalarROOTVGREAT;
    static const scalar SMALL = doubleScalarSMALL;
    static const scalar VSMALL = doubleScalarVSMALL;
    static const scalar ROOTVSMALL = doubleScalarROOTVSMALL;

    scalar readScalar(Istream& is);
}

#elif defined(WM_LDP)

// Define scalar as a long double

namespace Foam
{
    typedef longDoubleScalar scalar;

    static const scalar GREAT = longDoubleScalarGREAT;
    static const scalar VGREAT = longDoubleScalarVGREAT;
    static const scalar ROOTVGREAT = longDoubleScalarROOTVGREAT;
    static const scalar SMALL = longDoubleScalarSMALL;
    static const scalar VSMALL = longDoubleScalarVSMALL;
    static const scalar ROOTVSMALL = longDoubleScalarROOTVSMALL;

    scalar readScalar(Istream& is);
}

#else

#error "Undefined precision for scalar in scalar.H"

#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
