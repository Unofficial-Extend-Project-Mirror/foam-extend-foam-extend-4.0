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

Class
    faceSetAlgorithm

Description
    Class for convexSetAlgorithms pertaining to faces

Author
    Sandeep Menon
    University of Massachusetts Amherst
    All rights reserved

\*---------------------------------------------------------------------------*/

#include "faceSetAlgorithm.H"

#include "meshOps.H"
#include "triFace.H"
#include "polyMesh.H"
#include "triIntersection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
faceSetAlgorithm::faceSetAlgorithm
(
    const polyMesh& mesh,
    const pointField& newPoints,
    const UList<edge>& newEdges,
    const UList<face>& newFaces,
    const UList<cell>& newCells,
    const UList<label>& newOwner,
    const UList<label>& newNeighbour
)
:
    convexSetAlgorithm
    (
        mesh,
        newPoints,
        newEdges,
        newFaces,
        newCells,
        newOwner,
        newNeighbour
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void faceSetAlgorithm::computeNormFactor(const label index) const
{
    // Clear existing fields
    parents_.clear();

    centres_.clear();
    weights_.clear();

    // Compute refNorm and normFactor
    refNorm_ = newFaces_[index].normal(newPoints_);
    normFactor_ = mag(refNorm_);

    // Normalize for later use
    refNorm_ /= normFactor_ + VSMALL;

    // Compute a bounding box around the face
    box_ = boundBox(newFaces_[index].points(newPoints_), false);

    vector minToXb = (box_.min() - box_.midpoint());
    vector maxToXb = (box_.max() - box_.midpoint());

    // Scale it by a bit
    box_.min() += (1.5 * minToXb);
    box_.max() += (1.5 * maxToXb);
}


// Check whether the bounding box contains the entity
bool faceSetAlgorithm::contains(const label index) const
{
    // Fetch old face
    const face& checkFace = mesh_.faces()[index];

    // Check if the bounding box contains any of the supplied points
    forAll(checkFace, pointI)
    {
        if (box_.contains(mesh_.points()[checkFace[pointI]]))
        {
            return true;
        }
    }

    return false;
}


// Compute intersection
bool faceSetAlgorithm::computeIntersection
(
    const label newIndex,
    const label oldIndex,
    const label offset,
    bool output
) const
{
    bool intersects = false;

    const pointField& newPoints = newPoints_;
    const pointField& oldPoints = mesh_.points();

    const face& newFace = newFaces_[newIndex];
    const face& oldFace = mesh_.faces()[oldIndex];

    // Compute the normal for the old face
    vector oldNorm = oldFace.normal(oldPoints);

    // Normalize
    oldNorm /= mag(oldNorm) + VSMALL;

    if ((oldNorm & refNorm_) < 0.0)
    {
        // Opposite face orientation. Skip it.
        return false;
    }

    // Check if decomposition is necessary
    if (oldFace.size() > 3 || newFace.size() > 3)
    {
        // Decompose new / old faces
        DynamicList<FixedList<point, 3> > clippingTris(15);
        DynamicList<FixedList<point, 3> > subjectTris(15);

        label ntOld = 0, ntNew = 0;

        if (oldFace.size() == 3)
        {
            // Add a new entry
            subjectTris.append(FixedList<point, 3>(vector::zero));

            subjectTris[ntOld][0] = oldPoints[oldFace[0]];
            subjectTris[ntOld][1] = oldPoints[oldFace[1]];
            subjectTris[ntOld][2] = oldPoints[oldFace[2]];

            ntOld++;
        }
        else
        {
            // Configure tris from oldFace
            vector ofCentre = oldFace.centre(oldPoints);

            forAll(oldFace, pI)
            {
                // Add a new entry
                subjectTris.append(FixedList<point, 3>(vector::zero));

                subjectTris[ntOld][0] = oldPoints[oldFace[pI]];
                subjectTris[ntOld][1] = oldPoints[oldFace.nextLabel(pI)];
                subjectTris[ntOld][2] = ofCentre;

                ntOld++;
            }
        }

        if (newFace.size() == 3)
        {
            // Add a new entry
            clippingTris.append(FixedList<point, 3>(vector::zero));

            clippingTris[ntNew][0] = newPoints[newFace[0]];
            clippingTris[ntNew][1] = newPoints[newFace[1]];
            clippingTris[ntNew][2] = newPoints[newFace[2]];

            ntNew++;
        }
        else
        {
            // Configure tris from newFace
            vector nfCentre = newFace.centre(newPoints);

            forAll(newFace, pI)
            {
                // Add a new entry
                clippingTris.append(FixedList<point, 3>(vector::zero));

                clippingTris[ntNew][0] = newPoints[newFace[pI]];
                clippingTris[ntNew][1] = newPoints[newFace.nextLabel(pI)];
                clippingTris[ntNew][2] = nfCentre;

                ntNew++;
            }
        }

        // Accumulate area / centroid over all intersections
        bool foundIntersect = false;

        scalar totalArea = 0.0;
        vector totalCentre = vector::zero;

        // Loop through all clipping tris
        forAll(clippingTris, i)
        {
            // Initialize the intersector
            triIntersection intersector(clippingTris[i]);

            // Test for intersection and evaluate
            // against all subject tris
            forAll(subjectTris, j)
            {
                intersects = intersector.evaluate(subjectTris[j]);

                if (intersects)
                {
                    scalar area = 0.0;
                    vector centre = vector::zero;

                    // Fetch area and centre
                    intersector.getAreaAndCentre(area, centre);

                    // Accumulate area / centroid
                    totalArea += area;
                    totalCentre += (area * centre);

                    foundIntersect = true;

                    if (output)
                    {
                        writeVTK
                        (
                            "polygonIntersect_"
                          + Foam::name(newIndex)
                          + '_'
                          + Foam::name(oldIndex)
                          + '_' + Foam::name(i)
                          + '_' + Foam::name(j),
                            intersector.getIntersection()
                        );
                    }
                }
            }
        }

        // Size-up the internal lists
        if (foundIntersect && !output)
        {
            // Normalize centre
            totalCentre /= totalArea + VSMALL;

            // Normalize and check if this is worth it
            if (mag(totalArea/normFactor_) > SMALL)
            {
                // Size-up the internal lists
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(totalArea, weights_);
                meshOps::sizeUpList(totalCentre, centres_);

                intersects = true;
            }
            else
            {
                intersects = false;
            }
        }
    }
    else
    {
        // Configure points for clipping triangle
        FixedList<point, 3> clippingTri(vector::zero);

        // Fill in points
        clippingTri[0] = newPoints[newFace[0]];
        clippingTri[1] = newPoints[newFace[1]];
        clippingTri[2] = newPoints[newFace[2]];

        // Configure points for subject triangle
        FixedList<point, 3> subjectTri(vector::zero);

        // Fill in points
        subjectTri[0] = oldPoints[oldFace[0]];
        subjectTri[1] = oldPoints[oldFace[1]];
        subjectTri[2] = oldPoints[oldFace[2]];

        // Initialize the intersector
        triIntersection intersector(clippingTri);

        // Test for intersection and evaluate
        intersects = intersector.evaluate(subjectTri);

        if (intersects)
        {
            scalar area;
            vector centre;

            // Fetch area and centre
            intersector.getAreaAndCentre(area, centre);

            // Normalize and check if this is worth it
            if (mag(area/normFactor_) > SMALL)
            {
                if (output)
                {
                    writeVTK
                    (
                        "triIntersectNew_"
                      + Foam::name(newIndex),
                        List<FixedList<point, 3> >(1, clippingTri)
                    );

                    writeVTK
                    (
                        "triIntersectOld_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        List<FixedList<point, 3> >(1, subjectTri)
                    );

                    writeVTK
                    (
                        "triIntersect_"
                      + Foam::name(newIndex)
                      + '_'
                      + Foam::name(oldIndex),
                        intersector.getIntersection()
                    );
                }

                // Size-up the internal lists
                meshOps::sizeUpList((oldIndex - offset), parents_);
                meshOps::sizeUpList(area, weights_);
                meshOps::sizeUpList(centre, centres_);
            }
            else
            {
                intersects = false;
            }
        }
    }

    return intersects;
}


//- Write out tris as a VTK
void faceSetAlgorithm::writeVTK
(
    const word& name,
    const List<FixedList<point, 3> >& triList
) const
{
    // Fill up all points
    label pI = 0;

    List<face> allTris(triList.size(), face(3));
    pointField allPoints(3 * triList.size());

    forAll(triList, triI)
    {
        allTris[triI][0] = pI;
        allPoints[pI++] = triList[triI][0];

        allTris[triI][1] = pI;
        allPoints[pI++] = triList[triI][1];

        allTris[triI][2] = pI;
        allPoints[pI++] = triList[triI][2];
    }

    // Write out in face-to-node addressing
    meshOps::writeVTK
    (
        mesh_,
        name,
        identity(triList.size()),
        2,
        allPoints,
        List<edge>(0),
        allTris,
        List<cell>(0),
        List<label>(0)
    );
}


} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
