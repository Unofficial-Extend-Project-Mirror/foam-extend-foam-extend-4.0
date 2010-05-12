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

#include "faMeshWriter.H"
#include "writeFuns.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::faMeshWriter::faMeshWriter
(
    const faMesh& aMesh,
    const bool binary,
    const fileName& fName
)
:
    aMesh_(aMesh),
    binary_(binary),
    os_(fName.c_str())
{
    // Write header
    writeFuns::writeHeader(os_, binary_, faMesh::typeName);
    os_ << "DATASET UNSTRUCTURED_GRID" << std::endl;

    // Write topology
    label nFaceVerts = 0;

    const faceList& pp = aMesh.faces();

    forAll(pp, faceI)
    {
        nFaceVerts += pp[faceI].size() + 1;
    }

    os_ << "POINTS " << aMesh.nPoints() << " float" << std::endl;

    DynamicList<floatScalar> ptField(3*aMesh.nPoints());

    writeFuns::insert(aMesh.points(), ptField);
    writeFuns::write(os_, binary_, ptField);

    os_ << "CELLS " << aMesh.nFaces() << ' ' << nFaceVerts
        << std::endl;

    DynamicList<label> vertLabels(nFaceVerts);
    DynamicList<label> faceTypes(nFaceVerts);

    forAll(pp, faceI)
    {
        const face& f = pp[faceI];

        const label fSize = f.size();
        vertLabels.append(fSize);

        writeFuns::insert(f, vertLabels);

        if (fSize == 3)
        {
            faceTypes.append(vtkTopo::VTK_TRIANGLE);
        }
        else if (fSize == 4)
        {
            faceTypes.append(vtkTopo::VTK_QUAD);
        }
        else
        {
            faceTypes.append(vtkTopo::VTK_POLYGON);
        }
    }

    writeFuns::write(os_, binary_, vertLabels);

    os_ << "CELL_TYPES " << aMesh.nFaces() << std::endl;

    writeFuns::write(os_, binary_, faceTypes);
}


// ************************************************************************* //
