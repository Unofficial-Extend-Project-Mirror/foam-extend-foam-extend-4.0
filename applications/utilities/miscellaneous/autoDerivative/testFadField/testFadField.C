/*---------------------------------------------------------------------------*\
 =========                   |
 \\      /   F ield          | OpenFOAM: The Open Source CFD Toolbox
  \\    /    O peration      |
   \\  /     A nd            | Copyright held by original author
    \\/      M anipulation   |
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

Application
    testFadField

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "fadOneFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    fadScalarField tf(3, fadScalar(3.3));

    Info << "name: " << fadScalar::typeName << endl;

    Info << "tf value: " << FadOneValue(tf) << endl;
    Info << "tf deriv: " << FadOneDeriv(tf, 0) << endl;

    FadOneSetDeriv(tf, 0, Field<double>(3, 1));

    fadScalar ten(10.0);
    fadScalar bs = tf[0] - ten;
    fadScalar cs = tf[0] - 10;

    fadScalarField a = sqr(tf);
    fadScalarField b = tf - ten;
    fadScalarField c = ten - tf;
    fadScalarField d = sqr(tf - ten);
    fadScalarField e =
        (sqr(tf - bs + fadScalar(10)) + 22.3*sqr(cs - bs + fadScalar(10)));
    fadScalarField f = d/e;
//     fadScalarField g = Foam::atan(f);    // HJ, not prepared
//     fadScalarField h = Foam::log(f);     // HJ, not prepared


    Info << "new tf: " << tf << endl;

    Field<Vector<fadScalar> > tfVector
    (
        3,
        Vector<fadScalar>(fadScalar(1), fadScalar(2), fadScalar(3))
    );

    Info << "tfVector: " << tfVector << endl;

    return 0;
}


// ************************************************************************* //
