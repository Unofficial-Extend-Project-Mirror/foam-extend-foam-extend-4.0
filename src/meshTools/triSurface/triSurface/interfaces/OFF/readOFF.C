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
    Geomview OFF polyList format. Does triangulation.

\*---------------------------------------------------------------------------*/

#include "triSurface.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::readOFF(const fileName& OFFfileName)
{
    IFstream OFFfile(OFFfileName);

    if (!OFFfile.good())
    {
        FatalErrorIn("triSurface::readOFF(const fileName&)")
            << "Cannot read file " << OFFfileName
            << exit(FatalError);
    }

    // Read header
    string hdr = getLineNoComment(OFFfile);
    if (hdr != "OFF")
    {
        FatalErrorIn("triSurface::readOFF(const fileName&)")
            << "OFF file " << OFFfileName
            << " does not start with 'OFF'"
            << exit(FatalError);
    }


    label nPoints, nEdges, nElems;

    string line = getLineNoComment(OFFfile);
    {
        IStringStream lineStream(line);
        lineStream >> nPoints >> nElems >> nEdges;
    }

    // Read points
    pointField points(nPoints);

    forAll(points, pointi)
    {
        scalar x, y, z;
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);
            lineStream >> x >> y >> z;
        }
        points[pointi] = point(x, y, z);
    }

    // Read faces & triangulate them,
    DynamicList<labelledTri> tris(nElems);

    for (label faceI = 0; faceI < nElems; faceI++)
    {
        line = getLineNoComment(OFFfile);
        {
            IStringStream lineStream(line);

            label nVerts;
            lineStream >> nVerts;

            face f(nVerts);

            forAll(f, fp)
            {
                lineStream >> f[fp];
            }

            // Triangulate.
            if (nVerts == 3)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
            }
            else if (nVerts == 4)
            {
                tris.append(labelledTri(f[0], f[1], f[2], 0));
                tris.append(labelledTri(f[2], f[3], f[0], 0));
            }
            else
            {
                faceList triFaces(f.nTriangles(points));

                label nTri = 0;

                f.triangles(points, nTri, triFaces);

                forAll(triFaces, triFaceI)
                {
                    const face& f = triFaces[triFaceI];

                    tris.append(labelledTri(f[0], f[1], f[2], 0));
                }
            }
        }
    }

    tris.shrink();

    *this = triSurface(tris, points);

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
