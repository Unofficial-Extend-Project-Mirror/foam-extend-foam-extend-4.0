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
    CreateTurbulenceFields

Description
    Creates a full setturbulence fields.
    - Currently does not output nut and nuTilda

Source files:
    createFields.H

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "incompressible/singlePhaseTransportModel/singlePhaseTransportModel.H"
#include "incompressible/RASModel/RASModel.H"


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

#   include "createFields.H"

    for (label i=startTime; i<endTime; i++)
    {
        runTime.setTime(Times[i], i);

        Info<< "Time = " << runTime.timeName() << endl;

        // Cache the turbulence fields

        Info<< "\nRetrieving field k from turbulence model" << endl;
        const volScalarField k = RASModel->k();

        Info<< "\nRetrieving field epsilon from turbulence model" << endl;
        const volScalarField epsilon = RASModel->epsilon();

        Info<< "\nRetrieving field R from turbulence model" << endl;
        const volSymmTensorField R = RASModel->R();

        // Check availability of tubulence fields

        if (!IOobject("k", runTime.timeName(), mesh).headerOk())
        {
            Info<< "\nWriting turbulence field k" << endl;
            k.write();
        }
        else
        {
            Info<< "\nTurbulence k field already exists" << endl;
        }

        if (!IOobject("epsilon", runTime.timeName(), mesh).headerOk())
        {
            Info<< "\nWriting turbulence field epsilon" << endl;
            epsilon.write();
        }
        else
        {
            Info<< "\nTurbulence epsilon field already exists" << endl;
        }

        if (!IOobject("R", runTime.timeName(), mesh).headerOk())
        {
            Info<< "\nWriting turbulence field R" << endl;
            R.write();
        }
        else
        {
            Info<< "\nTurbulence R field already exists" << endl;
        }

        if (!IOobject("omega", runTime.timeName(), mesh).headerOk())
        {
            const scalar Cmu = 0.09;

            Info<< "creating omega" << endl;
            volScalarField omega
            (
                IOobject
                (
                    "omega",
                    runTime.timeName(),
                    mesh
                ),
                epsilon/(Cmu*k),
                epsilon.boundaryField().types()
            );
            Info<< "\nWriting turbulence field omega" << endl;
            omega.write();
        }
        else
        {
            Info<< "\nTurbulence omega field already exists" << endl;
        }
    }

    Info<< "\nEnd\n" << endl;

    return(0);
}


// ************************************************************************* //

