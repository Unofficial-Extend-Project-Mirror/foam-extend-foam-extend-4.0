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

\*---------------------------------------------------------------------------*/

#include "immersedBoundaryFvPatch.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::immersedBoundaryFvPatch::makeTriAddressing() const
{
    if (debug)
    {
        InfoIn("void immersedBoundaryFvPatch::makeTriAddressing() const")
            << "creating tri addressing for immersed boundary " << name()
            << endl;
    }

    // It is an error to attempt to recalculate
    // if the pointer is already set
    if (cellsToTriAddrPtr_ || cellsToTriWeightsPtr_)
    {
        FatalErrorIn("immersedBoundaryFvPatch::makeTriAddressing() const")
            << "tri addressing already exist"
            << "for immersed boundary" << name()
            << abort(FatalError);
    }

    // Get reference to tri patch and hit faces
    const triSurface& triPatch = ibPolyPatch_.ibMesh();
    const vectorField& triCentres = triPatch.faceCentres();

    const labelList& hf = hitFaces();
    const vectorField& ibp = ibPoints();

    // Create a markup field and mark all tris containing an ib point with its
    // index
    labelList hitTris(triPatch.size(), -1);

    forAll (hf, hfI)
    {
        hitTris[hf[hfI]] = hfI;
    }

    // Allocate storage
    cellsToTriAddrPtr_ = new labelListList(triPatch.size());
    labelListList& addr = *cellsToTriAddrPtr_;

    cellsToTriWeightsPtr_ = new scalarListList(triPatch.size());
    scalarListList& w = *cellsToTriWeightsPtr_;

    // Algorithm:
    // For each tri face, check if it contains an IB point
    // - if so, set the addressing to the index of IB point and weight to 1
    // - if not, search the neighbouring faces of the visited faces until
    //   at least 3 IB points are found, or the neighbourhood is exhausted.
    //   When a sufficient number of points is found, calculate the weights
    //   using inverse distance weighting

    // Get addressing from the triangular patch
    const labelListList& pf = triPatch.pointFaces();

    forAll (triPatch, triI)
    {
        if (hitTris[triI] > -1)
        {
            // Triangle contains IB point
            addr[triI].setSize(1);
            w[triI].setSize(1);

            addr[triI] = hitTris[triI];
            w[triI] = 1;
        }
        else
        {
            // No direct hit.  Start a neighbourhood search

            // Record already visited faces
            labelHashSet visited;

            // Collect new faces to visit
            SLList<label> nextToVisit;

            // Collect IB points for interpolation
            labelHashSet ibPointsToUse;

            // Initialise with the original tri
            nextToVisit.insert(triI);

            do
            {
                const label curTri = nextToVisit.removeHead();

                // Discard tri if already visited
                if (visited[curTri]) continue;

                visited.insert(curTri);

                const triFace& curTriPoints = triPatch[curTri];

                // For all current points of face, pick up neighbouring faces
                forAll (curTriPoints, tpI)
                {
                    const labelList curNbrs = pf[curTriPoints[tpI]];

                    forAll (curNbrs, nbrI)
                    {
                        if (!visited.found(curNbrs[nbrI]))
                        {
                            // Found a face which is not visited.  Add it to
                            // the list of faces to visit
                            nextToVisit.append(curNbrs[nbrI]);

                            if (hitTris[curNbrs[nbrI]] > -1)
                            {
                                // Found a neighbour with a hit: use this
                                // IB point
                                ibPointsToUse.insert(hitTris[curNbrs[nbrI]]);
                            }
                        }
                    }
                }
            } while
            (
                ibPointsToUse.size() < 3
             && !nextToVisit.empty()
            );

            // Found neighbourhood: collect addressing and weights
            addr[triI] = ibPointsToUse.toc();
            w[triI].setSize(addr[triI].size());

            labelList& curAddr = addr[triI];
            scalarList& curW = w[triI];

            vector curTriCentre = triCentres[triI];

            scalar sumW = 0;

            forAll (curAddr, ibI)
            {
                curW[ibI] = 1/mag(curTriCentre - ibp[curAddr[ibI]]);
                sumW += curW[ibI];
            }

            // Divide weights by sum distance
            forAll (curW, ibI)
            {
                curW[ibI] /= sumW;
            }
        }
    }
}


const Foam::labelListList&
Foam::immersedBoundaryFvPatch::cellsToTriAddr() const
{
    if (!cellsToTriAddrPtr_)
    {
        makeTriAddressing();
    }

    return *cellsToTriAddrPtr_;
}


const Foam::scalarListList&
Foam::immersedBoundaryFvPatch::cellsToTriWeights() const
{
    if (!cellsToTriWeightsPtr_)
    {
        makeTriAddressing();
    }

    return *cellsToTriWeightsPtr_;
}


// ************************************************************************* //
