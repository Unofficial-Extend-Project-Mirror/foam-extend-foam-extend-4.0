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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool Foam::moleculeCloud::testPointFaceDistance
(
    const label p,
    const label faceNo
) const
{
    const vector& pointPosition(mesh_.points()[p]);

    return testPointFaceDistance(pointPosition, faceNo);
}

bool Foam::moleculeCloud::testPointFaceDistance
(
    const label p,
    const referredCell& refCell
) const
{
    const vector& pointPosition(mesh_.points()[p]);

    forAll (refCell.faces(), rCF)
    {
        if
        (
            testPointFaceDistance
            (
                pointPosition,
                refCell.faces()[rCF],
                refCell.vertexPositions(),
                refCell.faceCentres()[rCF],
                refCell.faceAreas()[rCF]
            )
        )
        {
            return true;
        }
    }

    return false;
}

bool Foam::moleculeCloud::testPointFaceDistance
(
    const vectorList& pointsToTest,
    const label faceNo
) const
{
    forAll(pointsToTest, pTT)
    {
        const vector& p(pointsToTest[pTT]);

        // if any point in the list is in range of the face
        // then the rest do not need to be tested and
        // true can be returned

        if (testPointFaceDistance(p, faceNo))
        {
            return true;
        }
    }

    return false;
}

bool Foam::moleculeCloud::testPointFaceDistance
(
    const vector& p,
    const label faceNo
) const
{
    const face& faceToTest(mesh_.faces()[faceNo]);

    const vector& faceC(mesh_.faceCentres()[faceNo]);

    const vector& faceA(mesh_.faceAreas()[faceNo]);

    const vectorList& points(mesh_.points());

    return testPointFaceDistance
    (
        p,
        faceToTest,
        points,
        faceC,
        faceA
    );
}

bool Foam::moleculeCloud::testPointFaceDistance
(
    const vector& p,
    const labelList& faceToTest,
    const vectorList& points,
    const vector& faceC,
    const vector& faceA
) const
{
    vector faceN(faceA/mag(faceA));

    scalar perpDist((p - faceC) & faceN);

    if (mag(perpDist) > pairPotentials_.rCutMax())
    {
        return false;
    }

    vector pointOnPlane = (p - faceN * perpDist);

    if (magSqr(faceC - pointOnPlane) < pairPotentials_.rCutMaxSqr()*1e-8)
    {
        // If pointOnPlane is very close to the face centre
        // then defining the local axes will be inaccurate
        // and it is very likely that pointOnPlane will be
        // inside the face, so return true if the points
        // are in range to be safe

        return (magSqr(pointOnPlane - p) <= pairPotentials_.rCutMaxSqr());
    }

    vector xAxis = (faceC - pointOnPlane)/mag(faceC - pointOnPlane);

    vector yAxis =
        ((faceC - pointOnPlane) ^ faceN)
       /mag((faceC - pointOnPlane) ^ faceN);

    List<vector2D> local2DVertices(faceToTest.size());

    forAll(faceToTest, fTT)
    {
        const vector& V(points[faceToTest[fTT]]);

        if (magSqr(V-p) <= pairPotentials_.rCutMaxSqr())
        {
            return true;
        }

        local2DVertices[fTT] = vector2D
        (
            ((V - pointOnPlane) & xAxis),
            ((V - pointOnPlane) & yAxis)
        );
    }

    scalar localFaceCx((faceC - pointOnPlane) & xAxis);

    scalar la_valid = -1;

    forAll(local2DVertices, fV)
    {
        const vector2D& va(local2DVertices[fV]);

        const vector2D& vb
        (
            local2DVertices[(fV + 1) % local2DVertices.size()]
        );

        if (mag(vb.y()-va.y()) > SMALL)
        {
            scalar la =
                (
                    va.x() - va.y()*((vb.x() - va.x())/(vb.y() - va.y()))
                )
               /localFaceCx;

            scalar lv = -va.y()/(vb.y() - va.y());


            if (la >= 0 && la <= 1 && lv >= 0 && lv <= 1)
            {
                la_valid = la;

                break;
            }
        }
    }

    if (la_valid < 0)
    {
        // perpendicular point inside face, nearest point is pointOnPlane;
        return (magSqr(pointOnPlane-p) <= pairPotentials_.rCutMaxSqr());
    }
    else
    {
        // perpendicular point outside face, nearest point is
        // on edge that generated la_valid;
        return
        (
            magSqr(pointOnPlane + la_valid*(faceC - pointOnPlane) - p)
         <= pairPotentials_.rCutMaxSqr()
        );
    }

    // if the algorithm hasn't returned anything by now then something has
    // gone wrong.

    FatalErrorIn("moleculeCloudTestPointAndFaceDistance.C") << nl
        << "point " << p << " to face " << faceToTest
        << " comparison did not find a nearest point"
        << " to be inside or outside face."
        << abort(FatalError);

    return false;
}


// ************************************************************************* //
