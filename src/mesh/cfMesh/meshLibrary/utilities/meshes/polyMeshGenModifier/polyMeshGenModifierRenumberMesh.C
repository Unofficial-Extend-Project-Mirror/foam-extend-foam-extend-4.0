/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Author: Franjo Juretic (franjo.juretic@c-fields.com)
     \\/     M anipulation  | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of cfMesh.

    cfMesh is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    cfMesh is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with cfMesh.  If not, see <http://www.gnu.org/licenses/>.

Description

\*---------------------------------------------------------------------------*/

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"
#include "polyMeshGenAddressing.H"
#include "SLList.H"

//#define DEBUG_ZIPUP

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::renumberMesh()
{
    Info << "Renumbering the mesh" << endl;

    labelList newOrder(mesh_.cells().size());

    if( true )
    {
        const VRWGraph& cellCells = mesh_.addressingData().cellCells();

        //- the business bit of the renumbering
        labelLongList nextCell;

        boolList visited(cellCells.size(), false);

        label currentCell;
        label cellInOrder = 0;

        //- loop over the cells
        forAll(visited, cellI)
        {
            //- find the first cell that has not been visited yet
            if( !visited[cellI] )
            {
                currentCell = cellI;

                //- use this cell as a start
                nextCell.append(currentCell);

                //- loop through the nextCell list.
                //- Add the first cell into the
                //- cell order if it has not already been
                //- visited and ask for its
                //- neighbours. If the neighbour in question
                //- has not been visited,
                //- add it to the end of the nextCell list
                while( nextCell.size() > 0 )
                {
                    currentCell = nextCell.removeLastElement();

                    if( !visited[currentCell] )
                    {
                        visited[currentCell] = true;

                        //- add into cellOrder
                        newOrder[cellInOrder] = currentCell;
                        ++cellInOrder;

                        //- find if the neighbours have been visited
                        forAllRow(cellCells, currentCell, nI)
                        {
                            const label nei = cellCells(currentCell, nI);

                            if( !visited[nei] )
                            {
                                //- not visited, add to the list
                                nextCell.append(nei);
                            }
                        }
                    }
                }
            }
        }
    }

    cellListPMG& oldCells = this->cellsAccess();
    const labelList& oldOwner = mesh_.owner();
    const labelList& oldNeighbour = mesh_.neighbour();

    cellList newCells(oldCells.size());

    //- The reverse order list gives the new cell label for every old cell
    labelLongList reverseOrder(newOrder.size());

    forAll(newOrder, cellI)
    {
        newCells[cellI].transfer(oldCells[newOrder[cellI]]);

        reverseOrder[newOrder[cellI]] = cellI;
    }

    //- Renumber the faces.
    //- Reverse face order gives the new face number for every old face
    labelLongList reverseFaceOrder(oldOwner.size(), 0);

    //- Mark the internal faces with -2 so that they are inserted first
    forAll(newCells, cellI)
    {
        const cell& c = newCells[cellI];

        forAll(c, faceI)
        {
            --reverseFaceOrder[c[faceI]];
        }
    }

    //- Order internal faces
    label nMarkedFaces = 0;

    forAll(newCells, cellI)
    {
        //- Note:
        //- Insertion cannot be done in one go as the faces need to be
        //- added into the list in the increasing order of neighbour
        //- cells.  Therefore, all neighbours will be detected first
        //- and then added in the correct order.

        const cell& c = newCells[cellI];

        //- Record the neighbour cell
        DynList<label, 24> neiCells(c.size(), -1);

        label nNeighbours(0);

        forAll(c, faceI)
        {
            if( reverseFaceOrder[c[faceI]] == -2 )
            {
                //- Face is internal and gets reordered
                if( cellI == reverseOrder[oldOwner[c[faceI]]] )
                {
                    neiCells[faceI] = reverseOrder[oldNeighbour[c[faceI]]];
                }
                else if( cellI == reverseOrder[oldNeighbour[c[faceI]]] )
                {
                    neiCells[faceI] = reverseOrder[oldOwner[c[faceI]]];
                }
                else
                {
                    Info << "Screwed up!!!" << endl;
                }

                ++nNeighbours;
            }
        }

        //- Add the faces in the increasing order of neighbours
        for(label i=0;i<nNeighbours;++i)
        {
            //- Find the lowest neighbour which is still valid
            label nextNei = -1;
            label minNei = oldCells.size();

            forAll(neiCells, ncI)
            {
                if( (neiCells[ncI] > -1) && (neiCells[ncI] < minNei) )
                {
                    nextNei = ncI;
                    minNei = neiCells[ncI];
                }
            }

            if( nextNei > -1 )
            {
                //- Face is internal and gets reordered
                reverseFaceOrder[c[nextNei]] = nMarkedFaces;

                //- Stop the neighbour from being used again
                neiCells[nextNei] = -1;

                ++nMarkedFaces;
            }
            else
            {
                FatalErrorIn
                (
                    "void polyMeshGenModifier::renumberedMesh() const"
                )   << "Error in internal face insertion"
                    << abort(FatalError);
            }
        }
    }

    //- Insert the boundary faces into reordering list
    forAll(reverseFaceOrder, faceI)
    {
        if( reverseFaceOrder[faceI] < 0 )
        {
            reverseFaceOrder[faceI] = nMarkedFaces;

            ++nMarkedFaces;
        }
    }

    //- Face order gives the old face label for every new face
    labelLongList faceOrder(reverseFaceOrder.size());

    forAll(faceOrder, faceI)
    {
        faceOrder[reverseFaceOrder[faceI]] = faceI;
    }

    //- Renumber the cells
    forAll(newCells, cellI)
    {
        cell& c = newCells[cellI];

        forAll(c, fI)
        {
            c[fI] = reverseFaceOrder[c[fI]];
        }
    }

    faceListPMG& oldFaces = this->facesAccess();
    faceList newFaces(oldFaces.size());

    forAll(newFaces, faceI)
    {
        newFaces[faceI].transfer(oldFaces[faceOrder[faceI]]);
    }

    //- Turn the face that need to be turned
    //- Only loop through internal faces
    forAll(oldNeighbour, faceI)
    {
        const label oldFaceI = faceOrder[faceI];
        if( oldNeighbour[oldFaceI] < 0 )
            continue;

        if
        (
            reverseOrder[oldNeighbour[oldFaceI]]
          < reverseOrder[oldOwner[oldFaceI]]
        )
        {
            newFaces[faceI] = newFaces[faceI].reverseFace();
        }
    }

    //- transfer faces and cells back to the original lists
    forAll(newCells, cellI)
        oldCells[cellI].transfer(newCells[cellI]);
    forAll(newFaces, faceI)
        oldFaces[faceI].transfer(newFaces[faceI]);

    mesh_.updateFaceSubsets(reverseFaceOrder);
    mesh_.updateCellSubsets(reverseOrder);
    this->clearOut();
    mesh_.clearOut();

    Info << "Finished renumbering the mesh" << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
