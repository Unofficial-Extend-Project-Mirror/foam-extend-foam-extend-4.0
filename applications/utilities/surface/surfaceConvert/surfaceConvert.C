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

Description
    Converts to and from Foam surface format. Optionally orders triangles
    by region.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fileName.H"
#include "triSurface.H"
#include "OFstream.H"
#include "OSspecific.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validOptions.insert("cleanup", "");
    argList::validOptions.insert("group", "");
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList args(argc, argv);

    fileName inFileName(args.additionalArgs()[0]);
    fileName outFileName(args.additionalArgs()[1]);

    if (outFileName == inFileName)
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite input file."
            << exit(FatalError);
    }

    Info << "Reading : " << inFileName << endl;
    triSurface surf(inFileName);

    Info<< "Read surface:" << endl;
    surf.writeStats(Info);
    Info<< endl;
    

    if (args.options().found("cleanup"))
    {
        Info << "Cleaning up surface" << endl;
        surf.cleanup(true);

        Info<< "After cleaning up surface:" << endl;
        surf.writeStats(Info);
        Info<< endl;
    }

    bool sortByRegion = args.options().found("group");

    if (sortByRegion)
    {
        Info << "Reordering faces into groups; one per region." << endl;
    }
    else
    {
        Info << "Maintaining face ordering" << endl;
    }

    Info << "Writing : " << outFileName << endl;
    surf.write(outFileName, sortByRegion);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
