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

\*---------------------------------------------------------------------------*/

#include "molConfig.H"
#include "fvCFD.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    Info<< nl << "Reading molecular configuration description dictionary"
        << endl;

    IOobject molConfigDescriptionIOobject
    (
        "molConfigDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        false
    );

    if (!molConfigDescriptionIOobject.headerOk())
    {
        FatalErrorIn(args.executable())
            << "Cannot find molConfig description file " << nl
            << args.caseName()/runTime.system()/"molConfig"/"molConfigDict"
            << nl << exit(FatalError);
    }

    IOdictionary molConfigDescription(molConfigDescriptionIOobject);


    // Create molCloud, registering object with mesh

    Info<< nl << "Creating molecular configuration" << endl;

    molConfig molecules(molConfigDescription, mesh);

    label totalMolecules = molecules.nMol();

    if (Pstream::parRun())
    {
        reduce(totalMolecules, sumOp<label>());
    }

    Info<< nl << "Total number of molecules added: " << totalMolecules
        << nl << endl;

    moleculeCloud molCloud
    (
        mesh,
        molecules.nMol(),
        molecules.id(),
        molecules.mass(),
        molecules.positions(),
        molecules.cells(),
        molecules.U(),
        molecules.A(),
        molecules.tethered(),
        molecules.tetherPositions()
    );

    IOdictionary idListDict
    (
        IOobject
        (
            "idList",
            mesh.time().constant(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        )
    );

    idListDict.add("idList", molecules.molIdList());

    IOstream::defaultPrecision(12);

    Info << nl << "Writing molecular configuration" << endl;

    if (!mesh.write())
    {
        FatalErrorIn(args.executable())
            << "Failed writing moleculeCloud."
            << nl << exit(FatalError);
    }

    Info<< nl << "ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info << nl << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
