/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2004-2007 Z. Tukovic and H. Jasak
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
    finiteAreaFoam

Description

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "faCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

#   include "createFaMesh.H"
#   include "createFaFields.H"
#   include "createVolFields.H"

    Info<< "Total mass of surfactant: " 
        << sum(Cs.internalField()*aMesh.S()) << endl;

    Info << "\nStarting time loop\n" << endl;

    for (runTime++; !runTime.end(); runTime++)
    {
        Info << "Time = " << runTime.value() << endl;

        faScalarMatrix CsEqn
        (
            fam::ddt(Cs)
          + fam::div(phis, Cs)
          - fam::laplacian(Ds, Cs)
        );

        CsEqn.solve();

        if (runTime.outputTime())
        {
            vsm.mapToVolume(Cs, Cvf.boundaryField());

            runTime.write();
        }

        Info<< "Total mass of surfactant: " 
            << sum(Cs.internalField()*aMesh.S()) << endl;

        Info << "ExecutionTime = "
            << scalar(runTime.elapsedCpuTime())
            << " s\n" << endl << endl;
    }

    return(0);
}

// ************************************************************************* //
