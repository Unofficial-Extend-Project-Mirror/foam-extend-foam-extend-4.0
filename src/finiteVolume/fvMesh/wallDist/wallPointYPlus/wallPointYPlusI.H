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

Description

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// Update this with w2 if w2 nearer to pt.
inline bool wallPointYPlus::update
(
    const point& pt,
    const wallPointYPlus& w2,
    const scalar tol
)
{
    scalar dist2 = magSqr(pt - w2.origin());

    scalar diff = distSqr() - dist2;

    if (diff < 0)
    {
        // already nearer to pt
        return false;
    }

    if ((diff < SMALL) || ((distSqr() > SMALL) && (diff/distSqr() < tol)))
    {
        // don't propagate small changes
        return false;
    }
    else
    {
        // only propagate if interesting (i.e. y+ < 100)
        scalar yPlus = Foam::sqrt(dist2)/w2.data();


        if (yPlus < yPlusCutOff)
        {
            // update with new values
            distSqr() = dist2;
            origin() = w2.origin();
            data() = w2.data();

            return true;
        }
        else
        {
            return false;
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
inline wallPointYPlus::wallPointYPlus()
:
    wallPointData<scalar>()
{
    // Important: value of yStar where meshWave does not come.
    data() = 1.0;
}


// Construct from components
inline wallPointYPlus::wallPointYPlus
(
    const point& origin,
    const scalar yStar,
    const scalar distSqr
)
:
    wallPointData<scalar>(origin, yStar, distSqr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Update this with w2 if w2 nearer to pt.
inline bool wallPointYPlus::updateCell
(
    const polyMesh& mesh,
    const label thisCellI,
    const label neighbourFaceI,
    const wallPointYPlus& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& cellCentres = mesh.primitiveMesh::cellCentres();

    return update
    (
        cellCentres[thisCellI],
        neighbourWallInfo,
        tol
    );
}


// Update this with w2 if w2 nearer to pt.
inline bool wallPointYPlus::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const label neighbourCellI,
    const wallPointYPlus& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFaceI],
        neighbourWallInfo,
        tol
    );
}


// Update this with w2 if w2 nearer to pt.
inline bool wallPointYPlus::updateFace
(
    const polyMesh& mesh,
    const label thisFaceI,
    const wallPointYPlus& neighbourWallInfo,
    const scalar tol
)
{
    const vectorField& faceCentres = mesh.faceCentres();

    return update
    (
        faceCentres[thisFaceI],
        neighbourWallInfo,
        tol
    );
}


// ************************************************************************* //

} // End namespace Foam

// ************************************************************************* //
