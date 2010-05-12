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
    solidStress

Description
    Calculates and writes the scalar fields of the six components of the stress 
    tensor sigma for each time for linear stress analysis calculations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

#   include "addTimeOptions.H"
#   include "setRootCase.H"

#   include "createTime.H"

    // Get times list
    instantList Times = runTime.times();

    // set startTime and endTime depending on -time and -latestTime options
#   include "checkTimeOptions.H"

    runTime.setTime(Times[startTime], startTime);

#   include "createMesh.H"
#   include "readMechanicalProperties.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        mesh.readUpdate();

        IOobject Uheader
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ
        );

        // Check U exists
        if (Uheader.headerOk())
        {
            mesh.readUpdate();

            Info<< "    Reading U" << endl;
            volVectorField U(Uheader, mesh);

            volTensorField gradU = fvc::grad(U);

        volSymmTensorField sigma =
            rho*(2.0*mu*symm(gradU) + lambda*I*tr(gradU));

            volScalarField sigmaEq
            (
                IOobject
                (
                    "sigmaEq",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sqrt((3.0/2.0)*magSqr(dev(sigma)))
            );

            Info<< "Max sigmaEq = " << max(sigmaEq).value()
                << endl;
            sigmaEq.write();

            volScalarField sigmaxx
            (
                IOobject
                (
                    "sigmaxx",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::XX)
            );
            sigmaxx.write();

            volScalarField sigmayy
            (
                IOobject
                (
                    "sigmayy",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::YY)
            );
            sigmayy.write();

            volScalarField sigmazz
            (
                IOobject
                (
                    "sigmazz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::ZZ)
            );
            sigmazz.write();

            volScalarField sigmaxy
            (
                IOobject
                (
                    "sigmaxy",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::XY)
            );
            sigmaxy.write();

            volScalarField sigmaxz
            (
                IOobject
                (
                    "sigmaxz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::XZ)
            );
            sigmaxz.write();

            volScalarField sigmayz
            (
                IOobject
                (
                    "sigmayz",
                    runTime.timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                sigma.component(symmTensor::YZ)
            );
            sigmayz.write(); 
       }
        else
        {
            Info<< "    No U" << endl;
        }

        Info<< endl;
    }

    Info<< "End" << endl;

    return(0);
}


// ************************************************************************* //
