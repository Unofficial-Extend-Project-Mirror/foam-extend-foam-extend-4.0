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

#include "triSurfaceCleanupDuplicates.H"
#include "triSurfModifier.H"
#include "meshOctree.H"
#include "demandDrivenData.H"

# ifdef USE_OMP
#include <omp.h>
# endif
#include <set>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool triSurfaceCleanupDuplicates::checkDuplicateTriangles()
{
    labelLongList newTriangleLabel(surf_.size(), -1);

    const VRWGraph& pointTriangles = surf_.pointFacets();

    //- check if there exist duplicate triangles
    label counter(0);

    forAll(surf_, triI)
    {
        if( newTriangleLabel[triI] != -1 )
            continue;

        newTriangleLabel[triI] = counter;
        ++counter;

        const labelledTri& tri = surf_[triI];

        forAll(pointTriangles[tri[0]], ptI)
        {
            const label triJ = pointTriangles(tri[0], ptI);

            if( triJ <= triI )
                continue;

            const labelledTri& otherTri = surf_[triJ];

            if( tri == otherTri )
                newTriangleLabel[triJ] = newTriangleLabel[triI];
        }
    }

    Info << "Found " << (newTriangleLabel.size()-counter)
        << " duplicate triangles" << endl;

    //- return if there exist no duplicate triangles
    if( counter == newTriangleLabel.size() )
        return false;

    Info << "Current number of triangles" << surf_.size() << endl;
    Info << "New number of triangles " << counter << endl;

    //- create new list of triangles and store it in the surface mesh
    LongList<labelledTri> newTriangles(counter);

    forAll(newTriangleLabel, triI)
    {
        newTriangles[newTriangleLabel[triI]] = surf_[triI];
    }

    updateTriangleLabels(newTriangleLabel);

    triSurfModifier(surf_).facetsAccess().transfer(newTriangles);
    surf_.updateFacetsSubsets(newTriangleLabel);

    return true;
}

bool triSurfaceCleanupDuplicates::mergeDuplicatePoints()
{
    pointField& pts = const_cast<pointField&>(surf_.points());
    labelLongList newPointLabel(surf_.nPoints());
    bool foundDuplicates(false);

    # ifdef USE_OMP
    # pragma omp parallel
    # endif
    {
        # ifdef USE_OMP
        # pragma omp for
        # endif
        forAll(newPointLabel, pI)
            newPointLabel[pI] = pI;

        # ifdef USE_OMP
        # pragma omp barrier
        # endif
        //- check if there exist any vertices closer
        //- than the prescribed tolerance
        # ifdef USE_OMP
        # pragma omp for schedule(dynamic, 20)
        # endif
        for(label leafI=0;leafI<octree_.numberOfLeaves();++leafI)
        {
            DynList<label> ct;
            octree_.containedTriangles(leafI, ct);

            std::set<label> points;

            forAll(ct, ctI)
            {
                const label triI = newTriangleLabel_[ct[ctI]];

                if( triI < 0 )
                    continue;

                const labelledTri& tri = surf_[triI];

                forAll(tri, i)
                    points.insert(tri[i]);
            }

            for
            (
                std::set<label>::const_iterator it=points.begin();
                it!=points.end();
            )
            {
                const label pointI = *it;

                for
                (
                    std::set<label>::const_iterator nIt=++it;
                    nIt!=points.end();
                    ++nIt
                )
                    if( magSqr(pts[pointI] - pts[*nIt]) < sqr(tolerance_) )
                    {
                        foundDuplicates = true;
                        # ifdef USE_OMP
                        # pragma omp critical
                        # endif
                        newPointLabel[*nIt] = pointI;
                    }
            }
        }
    }

    //- find if there exist no duplicate points
    if( !foundDuplicates )
        return false;

    //- remove vertices and update node labels
    label counter(0);
    forAll(pts, pI)
        if( newPointLabel[pI] == pI )
        {
            newPointLabel[pI] = counter;
            if( counter < pI )
                pts[counter] = pts[pI];
            ++counter;
        }
        else
        {
            const label origI = newPointLabel[pI];
            newPointLabel[pI] = newPointLabel[origI];
        }

    Info << "Found " << (pts.size() - counter) << "duplicate points" << endl;

    pts.setSize(counter);

    //- remove triangles containing duplicate points
    LongList<labelledTri> newTriangles(surf_.facets());
    labelLongList newTriangleLabel(surf_.size(), -1);

    counter = 0;
    forAll(surf_, triI)
    {
        const labelledTri& tri = surf_[triI];
        const labelledTri newTri
        (
            newPointLabel[tri[0]],
            newPointLabel[tri[1]],
            newPointLabel[tri[2]],
            tri.region()
        );

        bool store(true);
        for(label i=0;i<2;++i)
            for(label j=i+1;j<3;++j)
                if( newTri[i] == newTri[j] )
                {
                    store = false;
                    break;
                }

        if( store )
        {
            newTriangles[counter] = newTri;
            newTriangleLabel[triI] = counter;
            ++counter;
        }
    }

    newTriangles.setSize(counter);

    //- update the surface
    triSurfModifier(surf_).facetsAccess().transfer(newTriangles);
    surf_.updateFacetsSubsets(newTriangleLabel);

    return true;
}

void triSurfaceCleanupDuplicates::updateTriangleLabels
(
    const labelLongList& newTriangleLabel
)
{
    //- update addressing between the original triangles and the cleaned mesh
    forAll(newTriangleLabel_, triI)
    {
        const label origI = newTriangleLabel_[triI];

        if( origI >= 0 )
            newTriangleLabel_[triI] = newTriangleLabel[origI];
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
