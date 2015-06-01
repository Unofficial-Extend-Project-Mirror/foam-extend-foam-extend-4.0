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
    Finds feature edges and corners of a triangulated surface

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "IFstream.H"
#include "fileName.H"
#include "triSurf.H"
#include "demandDrivenData.H"

#include "triSurfaceCopyParts.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
using namespace Foam;

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.clear();
    argList::validArgs.append("input surface file");
    argList::validArgs.append("output surface file");
    argList::validArgs.append("patch/subset name");
    argList args(argc, argv);

    const fileName inFileName(args.args()[1]);
    const fileName outFileName(args.args()[2]);

    if( outFileName == inFileName )
    {
        FatalErrorIn(args.executable())
            << "Output file " << outFileName
            << " would overwrite the input file."
            << exit(FatalError);
    }

    wordList patches(1);
    patches[0] = args.args()[3];

    triSurf originalSurface(inFileName);

    triSurfaceCopyParts copyPart(originalSurface);

    triSurf copy;

    copyPart.copySurface(patches, copy);

    copy.writeSurface(outFileName);

    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
