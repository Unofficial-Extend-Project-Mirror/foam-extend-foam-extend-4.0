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

Application
    ptot

Description
    For each time: calculate the total pressure.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#   include "createMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

        Info<< "Time = " << runTime.timeName() << endl;

        IOobject pheader
        (
            "p",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );


        // Check p and U exist
        if (pheader.headerOk() && Uheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading p" << endl;
            volScalarField p(pheader, mesh);

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            Info<< "    Calculating ptot" << endl;
            if (p.dimensions() == dimensionSet(0, 2, -2, 0, 0))
            {
                volScalarField ptot
                (
                    IOobject
                    (
                        "ptot",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ
                    ),
                    p + 0.5*magSqr(U)
                );
                ptot.write();
            }
            else
            {
                IOobject rhoheader
                (
                    "rho",
                    runTime.timeName(),
                    mesh,
                    IOobject::MUST_READ
                );

                // Check rho exists
                if (rhoheader.headerOk())
                {
                    Info<< "    Reading rho" << endl;
                    volScalarField rho(rhoheader, mesh);

                    volScalarField ptot
                    (
                        IOobject
                        (
                            "ptot",
                            runTime.timeName(),
                            mesh,
                            IOobject::NO_READ
                        ),
                        p + 0.5*rho*magSqr(U)
                    );
                    ptot.write();
                }
                else
                {
                    Info<< "    No rho" << endl;
                }
            }
        }
        else
        {
            Info<< "    No p or U" << endl;
        }

        Info<< endl;
    }

    return 0;
}


// ************************************************************************* //
