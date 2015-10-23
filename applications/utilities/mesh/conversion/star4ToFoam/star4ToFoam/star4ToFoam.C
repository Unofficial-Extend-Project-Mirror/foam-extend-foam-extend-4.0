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
    star4ToFoam

Description
    Converts a Star-CD (v4) pro-STAR mesh into FOAM format.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "foamTime.H"
#include "starMeshReader.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
//  Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validArgs.append("pro-STAR file prefix");
    argList::validOptions.insert("noscale", "");

#   include "setRootCase.H"
#   include "createTime.H"

    // rescale from [mm] to [m] by default
    scalar scaleFactor = 0.001;
    if (args.options().found("noscale"))
    {
        scaleFactor = 1.0;
    }

    // common error - forgot a trailing '.'
    string prefix(args.args()[3]);
    if (prefix[prefix.length() - 1] == '.')
    {
        prefix.resize(prefix.length()-1);
    }

    starMeshReader meshRead(fileName(prefix), runTime, scaleFactor);

    // Set the precision of the points data to 10
    IOstream::defaultPrecision(10);

    Info << "Writing mesh" << endl;
    meshRead.writeMesh();

    Info<< nl << "End\n" << endl;

    return(0);
}

// ************************************************************************* //
