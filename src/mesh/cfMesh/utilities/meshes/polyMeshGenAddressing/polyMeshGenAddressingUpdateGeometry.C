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

\*---------------------------------------------------------------------------*/

#include "polyMeshGenAddressing.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void polyMeshGenAddressing::updateGeometry
(
    const boolList& changedFace
)
{
    const pointFieldPMG& p = mesh_.points();
    const faceListPMG& faces = mesh_.faces();

    //- update face centres and face areas
    if( faceCentresPtr_ && faceAreasPtr_ )
    {
        vectorField& fCtrs = *faceCentresPtr_;
        vectorField& fAreas = *faceAreasPtr_;

        # ifdef USE_OMP
        # pragma omp parallel for if( faces.size() > 100 ) \
        schedule(dynamic, 10)
        # endif
        forAll(faces, faceI)
            if( changedFace[faceI] )
            {
                const face& f = faces[faceI];
                const label nPoints = f.size();

                // If the face is a triangle, do a direct calculation for
                // efficiency and to avoid round-off error-related problems
                if (nPoints == 3)
                {
                    fCtrs[faceI] = (1.0/3.0)*(p[f[0]] + p[f[1]] + p[f[2]]);
                    fAreas[faceI] =
                        0.5*((p[f[1]] - p[f[0]])^(p[f[2]] - p[f[0]]));
                }
                else
                {
                    vector sumN = vector::zero;
                    scalar sumA = 0.0;
                    vector sumAc = vector::zero;

                    point fCentre = p[f[0]];
                    for(label pI=1;pI<nPoints;++pI)
                    {
                        fCentre += p[f[pI]];
                    }

                    fCentre /= nPoints;

                    for(label pI=0;pI<nPoints;++pI)
                    {
                        const point& nextPoint = p[f.nextLabel(pI)];

                        vector c = p[f[pI]] + nextPoint + fCentre;
                        vector n = (nextPoint - p[f[pI]])^(fCentre - p[f[pI]]);
                        scalar a = mag(n);

                        sumN += n;
                        sumA += a;
                        sumAc += a*c;
                    }

                    fCtrs[faceI] = (1.0/3.0)*sumAc/(sumA + VSMALL);
                    fAreas[faceI] = 0.5*sumN;
                }
            }
    }

    //- update cell centres and cell volumes
    if( cellCentresPtr_ && cellVolumesPtr_ && faceCentresPtr_ && faceAreasPtr_ )
    {
        const vectorField& fCtrs = *faceCentresPtr_;
        const vectorField& fAreas = *faceAreasPtr_;
        vectorField& cellCtrs = *cellCentresPtr_;
        scalarField& cellVols = *cellVolumesPtr_;

        const labelList& own = mesh_.owner();
        const cellListPMG& cells = mesh_.cells();

        # ifdef USE_OMP
        # pragma omp parallel for if( cells.size() > 100 ) \
        schedule(dynamic, 10)
        # endif
        forAll(cells, cellI)
        {
            const cell& c = cells[cellI];

            bool update(false);
            forAll(c, fI)
                if( changedFace[c[fI]] )
                {
                    update = true;
                    break;
                }

            if( update )
            {
                cellCtrs[cellI] = vector::zero;
                cellVols[cellI] = 0.0;

                //- estimate position of cell centre
                vector cEst(vector::zero);
                forAll(c, fI)
                    cEst += fCtrs[c[fI]];
                cEst /= c.size();

                forAll(c, fI)
                    if( own[c[fI]] == cellI )
                    {
                        // Calculate 3*face-pyramid volume
                        const scalar pyr3Vol =
                            max
                            (
                                fAreas[c[fI]] &
                                (
                                    fCtrs[c[fI]] -
                                    cEst
                                ),
                                VSMALL
                            );

                        // Calculate face-pyramid centre
                        const vector pc =
                            (3.0/4.0)*fCtrs[c[fI]] + (1.0/4.0)*cEst;

                        // Accumulate volume-weighted face-pyramid centre
                        cellCtrs[cellI] += pyr3Vol*pc;

                        // Accumulate face-pyramid volume
                        cellVols[cellI] += pyr3Vol;
                    }
                    else
                    {
                        // Calculate 3*face-pyramid volume
                        const scalar pyr3Vol =
                            max
                            (
                                fAreas[c[fI]] &
                                (
                                    cEst - fCtrs[c[fI]]
                                ),
                                VSMALL
                            );

                        // Calculate face-pyramid centre
                        const vector pc =
                            (3.0/4.0)*fCtrs[c[fI]] + (1.0/4.0)*cEst;

                        // Accumulate volume-weighted face-pyramid centre
                        cellCtrs[cellI] += pyr3Vol*pc;

                        // Accumulate face-pyramid volume
                        cellVols[cellI] += pyr3Vol;
                    }

                cellCtrs[cellI] /= cellVols[cellI];
                cellVols[cellI] /= 3.0;
            }
        }
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
