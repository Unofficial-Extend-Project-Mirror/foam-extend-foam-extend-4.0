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
    foamMeshToElmer

Description
    Reads an foam mesh and writes a CSC/Elmer .nodes/.elements/.boundary format.

Usage
    - foamMeshToElmer [OPTION] \n
    Reads an foam mesh and writes a CSC/Elmer .nodes/.elements/.boundary format.

    @param -scale \<factor\>\n
    Specify an alternative geometry scaling factor.
    The default is @b 1.

Note
    The cellTable information available in the files
    @c constant/cellTable and @c constant/polyMesh/cellTableId
    will be used if available. Otherwise the cellZones are used when
    creating the cellTable information.

See Also
    Foam::cellTable, Foam::meshWriter and Foam::meshWriters::Elmer

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "timeSelector.H"
#include "objectRegistry.H"
#include "foamTime.H"
#include "polyMesh.H"
#include "ElmerMeshWriter.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    argList::noParallel();
    timeSelector::addOptions();

    argList::validOptions.insert("scale", "factor");
    argList::validOptions.insert("exclude", "pattern");

#   include "setRootCase.H"
#   include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);


    // default: do not rescale
    scalar scaleFactor = 1;
    if (args.optionReadIfPresent("scale", scaleFactor))
    {
        if (scaleFactor <= 0)
        {
            scaleFactor = 1;
        }
    }

    wordRe excludePattern;
    args.optionReadIfPresent("exclude", excludePattern);
    Info<<"Exclude pattern: \"" << excludePattern << "\"." <<endl;
    excludePattern.compile(wordRe::DETECT_NOCASE);

#   include "createPolyMesh.H"


    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);

#       include "getTimeIndex.H"

        polyMesh::readUpdateState state = mesh.readUpdate();

        if (!timeI || state != polyMesh::UNCHANGED)
        {
            meshWriters::Elmer writer(mesh, excludePattern, scaleFactor);
            fileName dirName("elmerMesh");
            if (state != polyMesh::UNCHANGED)
            {
                dirName += '_' + runTime.timeName();
            }

            if (!writer.write(dirName))
            {
                break; // Conversion failed
            }
        }

        Info << nl << endl;
    }

    Info << "End\n" << endl;

    return 0;
}

// ************************************************************************* //
