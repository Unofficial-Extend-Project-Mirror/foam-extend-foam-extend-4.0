/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

#include "triSurface.H"
#include "mergePoints.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::stitchTriangles
(
    const pointField& rawPoints,
    const scalar tol,
    bool verbose
)
{
    // Merge points
    labelList pointMap;
    pointField newPoints;
    bool hasMerged = mergePoints(rawPoints, tol, verbose, pointMap, newPoints);

    if (hasMerged)
    {
        if (verbose)
        {
            Pout<< "stitchTriangles : Merged from " << rawPoints.size()
                << " points down to " << newPoints.size() << endl;
        }

        // Reset the triangle point labels to the unique points array
        label newTriangleI = 0;
        forAll(*this, i)
        {
            const labelledTri& tri = operator[](i);
            labelledTri newTri
            (
                pointMap[tri[0]],
                pointMap[tri[1]],
                pointMap[tri[2]],
                tri.region()
            );

            if
            (
                (newTri[0] != newTri[1])
             && (newTri[0] != newTri[2])
             && (newTri[1] != newTri[2])
            )
            {
                operator[](newTriangleI++) = newTri;
            }
            else if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removing triangle " << i
                    << " with non-unique vertices." << endl
                    << "    vertices   :" << newTri << endl
                    << "    coordinates:" << newTri.points(newPoints)
                    << endl;
            }
        }

        // If empty triangles are detected, remove them from the list
        if (newTriangleI != size())
        {
            if (verbose)
            {
                Pout<< "stitchTriangles : "
                    << "Removed " << size() - newTriangleI
                    << " triangles" << endl;
            }
            setSize(newTriangleI);

            // Possibly compact out any unused points (since used only
            // by triangles that have just been deleted)
            // Done in two passes to save memory (pointField)

            // 1. Detect only
            PackedBoolList pointIsUsed(newPoints.size());

            label nPoints = 0;

            forAll(*this, i)
            {
                const labelledTri& tri = operator[](i);

                forAll(tri, fp)
                {
                    label pointI = tri[fp];
                    if (pointIsUsed.set(pointI, 1))
                    {
                        nPoints++;
                    }
                }
            }

            if (nPoints != newPoints.size())
            {
                // 2. Compact.
                pointMap.setSize(newPoints.size());
                label newPointI = 0;
                forAll(pointIsUsed, pointI)
                {
                    if (pointIsUsed[pointI])
                    {
                        newPoints[newPointI] = newPoints[pointI];
                        pointMap[pointI] = newPointI++;
                    }
                }
                newPoints.setSize(newPointI);

                newTriangleI = 0;
                forAll(*this, i)
                {
                    const labelledTri& tri = operator[](i);
                    operator[](newTriangleI++) = labelledTri
                    (
                        pointMap[tri[0]],
                        pointMap[tri[1]],
                        pointMap[tri[2]],
                        tri.region()
                    );
                }
            }
        }

        // Set the coordinates to the merged ones
        storedPoints().transfer(newPoints);
    }

    return hasMerged;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
