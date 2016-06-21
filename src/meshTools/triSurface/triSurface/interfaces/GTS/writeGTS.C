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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void triSurface::writeGTS(const bool writeSorted, Ostream& os) const
{
    // Write header
    os  << "# GTS file" << endl
        << "# Regions:" << endl;

    labelList faceMap;

    surfacePatchList myPatches(calcPatches(faceMap));

    // Print patch names as comment
    forAll(myPatches, patchI)
    {
        os  << "#     " << patchI << "    "
            << myPatches[patchI].name() << endl;
    }
    os  << "#" << endl;


    const pointField& ps = points();

    os  << "# nPoints  nEdges  nTriangles" << endl
        << ps.size() << ' ' << nEdges() << ' ' << size() << endl;

    // Write vertex coords

    forAll(ps, pointi)
    {
        os  << ps[pointi].x() << ' '
            << ps[pointi].y() << ' '
            << ps[pointi].z() << endl;
    }

    // Write edges.
    // Note: edges are in local point labels so convert
    const edgeList& es = edges();
    const labelList& meshPts = meshPoints();

    forAll(es, edgei)
    {
        os  << meshPts[es[edgei].start()] + 1 << ' '
            << meshPts[es[edgei].end()] + 1 << endl;
    }

    // Write faces in terms of edges.
    const labelListList& faceEs = faceEdges();

    if (writeSorted)
    {
        label faceIndex = 0;
        forAll(myPatches, patchI)
        {
            for
            (
                label patchFaceI = 0;
                patchFaceI < myPatches[patchI].size();
                patchFaceI++
            )
            {
                const label faceI = faceMap[faceIndex++];

                const labelList& fEdges = faceEdges()[faceI];

                os  << fEdges[0] + 1 << ' '
                    << fEdges[1] + 1 << ' '
                    << fEdges[2] + 1 << ' '
                    << (*this)[faceI].region() << endl;
            }
        }
    }
    else
    {
        forAll(faceEs, faceI)
        {
            const labelList& fEdges = faceEdges()[faceI];

            os  << fEdges[0] + 1 << ' '
                << fEdges[1] + 1 << ' '
                << fEdges[2] + 1 << ' '
                << (*this)[faceI].region() << endl;
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
