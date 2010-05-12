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

\*----------------------------------------------------------------------------*/

#include "moleculeCloud.H"
#include "polyBoundaryMeshEntries.H"
#include "PstreamCombineReduceOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::moleculeCloud::buildCellInteractionLists()
{
// #       include "moleculeCloudCodeSnippets.H"

    List<DynamicList<label> > directInteractionList(mesh_.nCells());

#   include "moleculeCloudBuildDirectInteractionList.H"

    directInteractionList_.setSize(mesh_.nCells());

    forAll(directInteractionList, transDIL)
    {
        directInteractionList_[transDIL].transfer
        (
            directInteractionList[transDIL].shrink()
        );
    }

//     for sorting DILs

//     forAll(directInteractionList_, dIL)
//     {
//         sort(directInteractionList_[dIL]);
//     }

    cellOccupancy_.setSize(mesh_.nCells());

    DynamicList<referredCell> referredInteractionList;

#   include "moleculeCloudBuildReferredInteractionList.H"

//#   include "moleculeCloudBuildReferredInteractionList_old.H"

    referredInteractionList_.setSize
    (
        referredInteractionList.size()
    );

    forAll(referredInteractionList, rIL)
    {
        referredInteractionList_[rIL] = referredInteractionList[rIL];
    }

    referredInteractionList_.setRealCellsInRange();
}


// ************************************************************************* //
