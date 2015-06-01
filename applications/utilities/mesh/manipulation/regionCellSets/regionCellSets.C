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
    Create cellSets from disconnected regions.

\*---------------------------------------------------------------------------*/

#include "argList.H"
#include "fvCFD.H"
#include "regionSplit.H"
#include "cellSet.H"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
#   include "createMesh.H"

    regionSplit cellRegion(mesh);

    Info<< "Number of disconnected regions: " << cellRegion.nRegions() << endl;

    List<labelHashSet> lhs(cellRegion.nRegions());

    forAll (cellRegion, cellI)
    {
        lhs[cellRegion[cellI]].insert(cellI);
    }

    // Create sets
    forAll (lhs, setI)
    {
        word setName("cellRegion" + name(setI));

        Info<< "Creating set " << setName << " with " << lhs[setI].size()
            << " cells." << endl;

        cellSet
        (
            mesh,
            setName,
            lhs[setI]
        ).write();
    }

    Info << nl << "End" << endl;

    return 0;
}


// ************************************************************************* //
