/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "writePatchGeom.H"
#include "OFstream.H"
#include "floatScalar.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

void writePatchGeom
(
    const bool binary,
    const faceList& faces,
    const pointField& points,
    std::ofstream& pStream
)
{
    pStream << "POINTS " << points.size() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*points.size());

    writeFuns::insert(points, ptField);

    writeFuns::write(pStream, binary, ptField);


    label nFaceVerts = 0;

    forAll(faces, faceI)
    {
        nFaceVerts += faces[faceI].size() + 1;
    }
    pStream << "POLYGONS " << faces.size() << ' ' << nFaceVerts
        << std::endl;


    DynamicList<label> vertLabels(nFaceVerts);

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        vertLabels.append(f.size());

        writeFuns::insert(f, vertLabels);
    }
    writeFuns::write(pStream, binary, vertLabels);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
