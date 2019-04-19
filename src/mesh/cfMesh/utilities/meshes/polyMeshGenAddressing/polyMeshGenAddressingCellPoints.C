/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

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

void polyMeshGenAddressing::calcCellPoints() const
{
    if( cpPtr_ )
    {
        FatalErrorIn("polyMeshGenAddressing::calcCellPoints() const")
            << "cellPoints already calculated"
            << abort(FatalError);
    }
    else
    {
        const cellListPMG& cells = mesh_.cells();
        const faceListPMG& faces = mesh_.faces();

        //- create the storage
        cpPtr_ = new VRWGraph(cells.size());
        VRWGraph& cellPointsAddr = *cpPtr_;

        labelList nPoints(cells.size());

        # ifdef USE_OMP
        const label nThreads = 3 * omp_get_num_procs();
        # pragma omp parallel num_threads(nThreads) if( cells.size() > 10000 )
        # endif
        {
            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(nPoints, i)
                nPoints[i] = i;

            # ifdef USE_OMP
            # pragma omp for schedule(static)
            # endif
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label, 32> cPoints;
                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                        cPoints.appendIfNotIn(f[pI]);
                }

                nPoints[cellI] = cPoints.size();
            }

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp master
            # endif
            VRWGraphSMPModifier(cellPointsAddr).setSizeAndRowSize(nPoints);

            # ifdef USE_OMP
            # pragma omp barrier

            # pragma omp for schedule(static)
            # endif
            forAll(cells, cellI)
            {
                const cell& c = cells[cellI];

                DynList<label, 32> cPoints;
                forAll(c, fI)
                {
                    const face& f = faces[c[fI]];

                    forAll(f, pI)
                        cPoints.appendIfNotIn(f[pI]);
                }

                cellPointsAddr.setRow(cellI, cPoints);
            }
        }
    }
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const VRWGraph& polyMeshGenAddressing::cellPoints() const
{
    if( !cpPtr_ )
    {
        # ifdef USE_OMP
        if( omp_in_parallel() )
            FatalErrorIn
            (
                "const VRWGraph& polyMeshGenAddressing::cellPoints() const"
            ) << "Calculating addressing inside a parallel region."
                << " This is not thread safe" << exit(FatalError);
        # endif

        calcCellPoints();
    }

    return *cpPtr_;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
