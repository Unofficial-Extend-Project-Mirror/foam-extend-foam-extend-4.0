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

#include "vtk.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::vtk<Type>::writeGeometry
(
    const pointField& points,
    const faceList& faces,
    Ostream& os
) const
{
    // Write vertex coordinates

    os
        << "# vtk DataFile Version 2.0" << nl
        << "sampleSurface" << nl
        << "ASCII" << nl
        << "DATASET POLYDATA" << nl;

    os  << "POINTS " << points.size() << " float" << nl;

    forAll(points, pointI)
    {
        const point& pt = points[pointI];

        os  << float(pt.x()) << ' ' << float(pt.y()) << ' ' << float(pt.z())
            << nl;
    }
    os  << endl;


    // Write triangles

    label nFaceVerts = 0;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        nFaceVerts += f.size() + 1;
    }

    os  << "POLYGONS " << faces.size() << ' ' << nFaceVerts << nl;

    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        os << f.size();

        forAll(f, fp)
        {
            os << ' ' << f[fp];
        }
        os << nl;
    }
}


// Write scalarField in vtk format
template<class Type>
void Foam::vtk<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const scalarField& values,
    Ostream& os
) const
{
    // Write data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        os  << "CELL_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }

    os  << fieldName << " 1 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        os << float(values[elemI]);

        if (elemI > 0 && (elemI%10) == 0)
        {
            os << nl;
        }
        else
        {
            os << ' ';
        }
    }
    os << nl;
}


// Write vectorField in vtk format
template<class Type>
void Foam::vtk<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const vectorField& values,
    Ostream& os
) const
{
    // Write data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        os  << "CELL_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }

    os  << fieldName << " 3 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        const vector& v = values[elemI];

        os << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2]) << nl;
    }
}


// Write sphericalTensorField in vtk format
template<class Type>
void Foam::vtk<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const sphericalTensorField& values,
    Ostream& os
) const
{
    // Write data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        os  << "CELL_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }

    os  << fieldName << " 1 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        const sphericalTensor& v = values[elemI];

        os  << float(v[0])
            << nl;
    }
}


// Write symmTensorField in vtk format
template<class Type>
void Foam::vtk<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const symmTensorField& values,
    Ostream& os
) const
{
    // Write data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        os  << "CELL_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }

    os  << fieldName << " 6 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        const symmTensor& v = values[elemI];

        os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
            << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
            << nl;
    }
}


// Write tensorField in vtk format
template<class Type>
void Foam::vtk<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const tensorField& values,
    Ostream& os
) const
{
    // Write data
    if (values.size() == points.size())
    {
        os  << "POINT_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }
    else
    {
        os  << "CELL_DATA " << values.size()
            << nl
            << "FIELD attributes 1" << nl;
    }

    os  << fieldName << " 9 " << values.size() << " float" << nl;

    forAll(values, elemI)
    {
        const tensor& v = values[elemI];

        os  << float(v[0]) << ' ' << float(v[1]) << ' ' << float(v[2])
            << float(v[3]) << ' ' << float(v[4]) << ' ' << float(v[5])
            << float(v[6]) << ' ' << float(v[7]) << ' ' << float(v[8])
            << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
Foam::vtk<Type>::vtk()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::vtk<Type>::~vtk()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::vtk<Type>::write
(
    const fileName& samplePath,
    const fileName& timeDir,
    const fileName& surfaceName,
    const pointField& points,
    const faceList& faces,
    const fileName& fieldName,
    const Field<Type>& values,
    const bool verbose
) const
{
    fileName surfaceDir(samplePath/timeDir);

    if (!exists(surfaceDir))
    {
        mkDir(surfaceDir);
    }

    fileName planeFName(surfaceDir/fieldName + '_' + surfaceName + ".vtk");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << planeFName << endl;
    }

    OFstream vtkFile(planeFName);

    writeGeometry(points, faces, vtkFile);

    writeData(fieldName, points, values, vtkFile);
}


// ************************************************************************* //
