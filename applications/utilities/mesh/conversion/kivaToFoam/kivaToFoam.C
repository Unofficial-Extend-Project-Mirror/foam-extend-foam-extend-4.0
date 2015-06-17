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
    kivaToFoam

Description
    Converts a KIVA3v grid to FOAM format

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "IFstream.H"
#include "OFstream.H"
#include "cellShape.H"
#include "cellModeller.H"
#include "preservePatchTypes.H"
#include "emptyPolyPatch.H"
#include "wallPolyPatch.H"
#include "symmetryPolyPatch.H"
#include "wedgePolyPatch.H"
#include "cyclicPolyPatch.H"
#include "unitConversion.H"

using namespace Foam;

//- Supported KIVA versions
enum kivaVersions
{
    kiva3,
    kiva3v
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::noParallel();
    argList::validOptions.insert("file", "fileName");
    argList::validOptions.insert("version", "[kiva3|kiva3v]");
    argList::validOptions.insert("zHeadMin", "scalar");

#   include "setRootCase.H"
#   include "createTime.H"

    fileName kivaFileName("otape17");
    if (args.optionFound("file"))
    {
        kivaFileName = args.option("file");
    }

    kivaVersions kivaVersion = kiva3v;
    if (args.optionFound("version"))
    {
        word kivaVersionName = args.option("version");

        if (kivaVersionName == "kiva3")
        {
            kivaVersion = kiva3;
        }
        else if (kivaVersionName == "kiva3v")
        {
            kivaVersion = kiva3v;
        }
        else
        {
            FatalErrorIn("main(int argc, char *argv[])")
                << "KIVA file version " << kivaVersionName << " not supported"
                << exit(FatalError);

            args.printUsage();
            FatalError.exit(1);
        }
    }

    scalar zHeadMin = -GREAT;
    args.optionReadIfPresent("zHeadMin", zHeadMin);

#   include "readKivaGrid.H"

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
