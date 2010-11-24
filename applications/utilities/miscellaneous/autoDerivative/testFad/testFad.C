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
    testFad

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "FadOne.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    typedef FadOne<1> fadScalar;

    fadScalar x(1);
    x.deriv(0) = 1;

    fadScalar f1 = 2.0*Foam::pow(x, 2) + 3.0*x + 4.0;
    Info << "fadScalar of f1 is: "
         << f1 << " It should be: [9 7]" << endl;

    fadScalar y(1);
    y.deriv(0) = 1;

    fadScalar f2 = 3*Foam::pow(y, 2)*Foam::sin(Foam::pow(y, 2));
    Info << "fadScalar of f2 is: " << f2
         << " It should be: [2.5244129 8.2906397]" << endl;

    fadScalar z(4);
    z.deriv(0) = 1;

    fadScalar f3 = Foam::sqrt(Foam::pow(z, 3));
    Info << "fadScalar of f3 is: " << f3
         << " It should be: [8 3]" << endl;

    fadScalar f4 = Foam::atan(x);
    Info << "fadScalar of f4 is: " << f4
         << " It should be: [?, ?]" << endl;

    fadScalar f5 = Foam::log(x);
    Info << "fadScalar of f4 is: " << f5
         << " It should be: [?, ?]" << endl;

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
