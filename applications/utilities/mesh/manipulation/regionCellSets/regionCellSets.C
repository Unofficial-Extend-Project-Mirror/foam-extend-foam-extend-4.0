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
