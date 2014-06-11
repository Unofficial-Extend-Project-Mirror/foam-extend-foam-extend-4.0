/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | cfMesh: A library for mesh generation
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "VRWGraphSMPModifier.H"

# ifdef USE_OMP
#include <omp.h>
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void polyMeshGenAddressing::calcCellCells() const
{
    if( ccPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellCells() const")
            << "cellCells already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellListPMG& cells = mesh_.cells();

        const labelList& own = mesh_.owner();
        const labelList& nei = mesh_.neighbour();

        //- create the storage
        ccPtr_ = new VRWGraph();
        VRWGraph& cellCellAddr = *ccPtr_;

        labelList nNei(cells.size());

        # ifdef USE_OMP
        const label nThreads = 3 * omp_get_num_procs();
        # endif

        # ifdef USE_OMP
        # pragma omp parallel num_threads(nThreads)
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(nNei, i)
                nNei[i] = 0;

            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label> neiCells;

                forAll(c, fI)
                {
                    label neiCell = own[c[fI]];
                    if( (neiCell == cellI) && (nei[c[fI]] != -1) )
                        neiCell = nei[c[fI]];

                    if( neiCell != cellI )
                        neiCells.appendIfNotIn(neiCell);
                }

                nNei[cellI] = neiCells.size();
            }

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp master
            # endif
            VRWGraphSMPModifier(cellCellAddr).setSizeAndRowSize(nNei);

            # ifdef USE_OMP
            # pragma omp barrier

            //- fill the graph with data
            # pragma omp for schedule(static)
            # endif
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label> neiCells;

                forAll(c, fI)
                {
                    label neiCell = own[c[fI]];
                    if( (neiCell == cellI) && (nei[c[fI]] != -1) )
                        neiCell = nei[c[fI]];

                    if( neiCell != cellI )
                        neiCells.appendIfNotIn(neiCell);
                }

                cellCellAddr.setRow(cellI, neiCells);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::cellCells() const
{
    if( !ccPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::cellCells() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellCells();
    }

    return *ccPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
