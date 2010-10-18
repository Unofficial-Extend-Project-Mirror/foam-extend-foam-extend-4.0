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

\*---------------------------------------------------------------------------*/

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline scalar lengthScaleEstimator::ratioMin() const
{
    return ratioMin_;
}


inline scalar lengthScaleEstimator::ratioMax() const
{
    return ratioMax_;
}


inline scalar lengthScaleEstimator::growthFactor() const
{
    return growthFactor_;
}


//- Limit length scale for surface-edges
inline void lengthScaleEstimator::limitScale(scalar& scale) const
{
    // Need to account for multiplication by ratioMin / ratioMax as well
    scale = Foam::max(scale, minLengthScale_ / ratioMin_);
    scale = Foam::min(scale, maxLengthScale_ / ratioMax_);
}


//- Check if a particular patch is free-floating
//  (i.e., gets its length-scale from the interior)
inline bool lengthScaleEstimator::isFreePatch
(
    const label pIndex
) const
{
    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    if (freePatches_.found(boundary[pIndex].name()))
    {
        return true;
    }

    // Certain boundary types are also considered floating
    const word& pType = boundary[pIndex].type();

    if
    (
        (pType == "processor") ||
        (pType == "cyclic") ||
        (pType == "wedge") ||
        (pType == "empty") ||
        (pType == "symmetryPlane")
    )
    {
        return true;
    }

    return false;
}


//- Check if a particular patch is flagged
//  for proximity-based refinement
inline bool lengthScaleEstimator::isProximityPatch
(
    const label pIndex
) const
{
    if (proximityPatches_.found(mesh_.boundaryMesh()[pIndex].name()))
    {
        return true;
    }

    return false;
}


//- Check if a particular patch is flagged
//  for curvature-based refinement
inline bool lengthScaleEstimator::isCurvaturePatch
(
    const label pIndex
) const
{
    if (curvaturePatches_.found(mesh_.boundaryMesh()[pIndex].name()))
    {
        return true;
    }

    return false;
}


//- Check whether a particular point is too close
//  to a previous mesh slice location
inline bool lengthScaleEstimator::checkOldSlices
(
    const vector& gCentre
) const
{
    forAll(sliceBoxes_, boxI)
    {
        if (sliceBoxes_[boxI].contains(gCentre))
        {
            // Too close to another slice-point. Bail out.
            return true;
        }
    }

    return false;
}


//- Add a boundBox to the existing set of sliceBoxes
inline void lengthScaleEstimator::appendBox
(
    const boundBox& bBox
)
{
    label currentSize = sliceBoxes_.size();
    sliceBoxes_.setSize(currentSize + 1);
    sliceBoxes_[currentSize] = bBox;
}


//- Clear the list of sliceBoxes
inline void lengthScaleEstimator::clearBoxes()
{
    sliceBoxes_.clear();
}


//- Check whether a particular patch permits refinement
//  - Return true if not permitted
inline bool lengthScaleEstimator::checkRefinementPatch
(
    const label pIndex
) const
{
    if (pIndex < 0)
    {
        return false;
    }

    if (findIndex(noModPatchIDs_, pIndex) > -1)
    {
        return true;
    }

    // Refinement is allowed
    return false;
}


// Return the fixed length-scale value for a boundary face
inline scalar lengthScaleEstimator::fixedLengthScale
(
    const label fIndex,
    const label pIndex,
    bool usePolyMesh
) const
{
    scalar scale = 0.0;

    const polyBoundaryMesh& boundary = mesh_.boundaryMesh();

    // Check fixed length-scale patches
    // If the value is negative, average face length-scales.
    if (fixedPatches_.found(boundary[pIndex].name()))
    {
        scalar dictValue =
        (
            fixedPatches_[boundary[pIndex].name()][0].scalarToken()
        );

        if (dictValue > 0.0)
        {
            return dictValue;
        }
    }

    // Approximate a length-scale from face area
    if (usePolyMesh)
    {
        scale = Foam::sqrt(2.0 * mag(mesh_.faceAreas()[fIndex]));
    }

    return scale;
}


// Test an edge / face for proximity with other faces on proximity patches
// and return the scalar distance to an oppositely-oriented face.
inline bool lengthScaleEstimator::testProximity
(
    const vector& gCentre,
    const vector& gNormal,
    const scalar testStep,
    label& proxFace,
    scalar& proxDistance
) const
{
    // Reset input
    proxDistance = GREAT;

    DynamicList<label> posIndices(20);
    scalar minDeviation = -0.9;
    label nD = spatialRes_, binSize = proximityBins_.size();

    // Obtain min/max extents
    const point& bMin = proxBoundBox_.min();
    const point& bMax = proxBoundBox_.max();

    // Extend bounding-box dimensions a bit to avoid edge-effects.
    scalar ext = 0.02*(mag(bMax - bMin));

    // Define an inverse grid-cell size.
    scalar xL = nD/(bMax.x() - bMin.x() + ext);
    scalar yL = nD/(bMax.y() - bMin.y() + ext);
    scalar zL = nD/(bMax.z() - bMin.z() + ext);

    // Reset the proximity face
    proxFace = -1;

    // Now take multiple steps in both normal directions,
    // and add to the list of boxes to be checked.
    for (scalar dir = -1.0; dir < 2.0; dir += 2.0)
    {
        for (scalar step = 0.0; step < 5.0*testStep; step += testStep)
        {
            // Hash the point-location
            point p = (gCentre + (dir*step*gNormal)) - bMin;

            label i = label(mag(::floor(p.x()*xL)));
            label j = label(mag(::floor(p.y()*yL)));
            label k = label(mag(::floor(p.z()*zL)));

            label pos = mag(((k*nD*nD)+(j*nD)+i) % binSize);

            if (findIndex(posIndices, pos) == -1)
            {
                posIndices.append(pos);
            }
        }
    }

    // Obtain old-mesh face geometry for reference.
    const vectorField& faceAreas = mesh_.faceAreas();
    const vectorField& faceCentres = mesh_.faceCentres();

    forAll(posIndices, indexI)
    {
        const labelList& posBin = proximityBins_[posIndices[indexI]];

        forAll(posBin, faceI)
        {
            // Step 1: Measure the distance to the face.
            vector rFace = (faceCentres[posBin[faceI]] - gCentre);

            scalar distance = mag(rFace);

            // Step 2: Check if this face is oriented away from face / edge.
            const vector& fNormal = faceAreas[posBin[faceI]];

            scalar deviation = (gNormal & (fNormal/(mag(fNormal) + VSMALL)));

            if
            (
                (deviation < minDeviation) &&
                (distance < proxDistance)
            )
            {
                // Update statistics
                proxFace = posBin[faceI];
                proxDistance = distance;
                // minDeviation = deviation;
            }
        }
    }

    // Check for fall below threshold
    if ((proxFace != -1) && (proxDistance < sliceThreshold_))
    {
        return true;
    }

    // Does not fall below threshold
    return false;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
