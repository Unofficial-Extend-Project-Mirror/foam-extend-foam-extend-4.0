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

Application
    PODsolver

Description
    Proper orthogonal decomposition solver with run-time selectable
    physical model

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PODODE.H"
#include "ODESolver.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "setRootCase.H"

#   include "createTime.H"
#   include "createMesh.H"

    Info<< "Reading PODsolverDict\n" << endl;

    IOdictionary PODsolverDict
    (
        IOobject
        (
            "PODsolverDict",
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Create equation
    autoPtr<PODODE> equation
    (
        PODODE::New(mesh, PODsolverDict)
    );

    // Create solver
    autoPtr<ODESolver> solver =
        ODESolver::New
        (
            PODsolverDict.lookup("solver"),
            equation()
        );

    // Read required accuracy
    scalar eps = readScalar(PODsolverDict.lookup("eps"));

    // Interpolation coefficients
    scalarField& PODcoeffs = equation->coeffs();

    Info<< "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        solver->solve
        (
            runTime.value(),
            runTime.value() + runTime.deltaT().value(),
            eps,
            runTime.deltaT().value()
        );

        // Create fields and recalculate properties
        if (runTime.outputTime())
        {
            Info << "Writing fields" << endl;
            equation->write();
        }
    }

    Info<< "End\n" << endl;

    return(0);
}


// ************************************************************************* //
