/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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
    Set inital film thickness for droplet spreading case.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"
#   include "createFaMesh.H"
#   include "readTransportProperties.H"

    Info << "Reading field h" << endl;
    areaScalarField h
    (
        IOobject
        (
            "h",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        aMesh
    );

    scalar D = 0.000038;
    scalar A = 0.0014;
    scalar R = (sqr(A) + sqr(D))/(2*D);

    Info << "Spherical cap radius: " << R << " m" << endl;

    scalarField& hI = h.internalField();

    const vectorField& Cf = aMesh.areaCentres().internalField();

    scalarField a = sqrt
    (
        sqr(Cf.component(vector::X))
      + sqr(Cf.component(vector::Y))
    );

    hI = pos(A - a)*(sqrt(sqr(R) - sqr(a)) - (R - D))
       + neg(A - a)*h0.value();

    h.write();

    Info<< "\nEnd" << endl;

    return 0;
}


// ************************************************************************* //
