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

#include "rawSurfaceWriter.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"
#include "IOmanip.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::rawSurfaceWriter<Type>::writeGeometry
(
    const pointField& points,
    const label& pointI,
    Ostream& os
) const
{
    const point& pt = points[pointI];

    os << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ';
}


template<class Type>
void Foam::rawSurfaceWriter<Type>::writeGeometry
(
    const pointField& points,
    const faceList& faces,
    const label& faceI,
    Ostream& os
) const
{
    const point& ct = faces[faceI].centre(points);

    os << ct.x() << ' ' << ct.y() << ' ' << ct.z() << ' ';
}

// Write scalarField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const scalarField& values,
    Ostream& os
) const
{
    // header
    os << "# " << fieldName;

    if (values.size() == points.size())
    {
        os  << "  POINT_DATA " << values.size()
            << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size()
            << nl;
    }

    os  << "#  x  y  z  " << fieldName
        << endl;

    // Write data
    forAll(values, elemI)
    {
        if (values.size() == points.size())
        {
            writeGeometry(points, elemI, os);
        }
        else
        {
            writeGeometry(points, faces, elemI, os);
        }
        os << values[elemI] << endl;
    }
    os << nl;
}


// Write vectorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const vectorField& values,
    Ostream& os
) const
{
    // header
    os << "# " << fieldName;

    if (values.size() == points.size())
    {
        os  << "  POINT_DATA " << values.size()
            << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size()
            << nl;
    }

    os  << "#  x  y  z  "
        << fieldName << "_x  "
        << fieldName << "_y  "
        << fieldName << "_z  "
        << endl;

    // Write data
    forAll(values, elemI)
    {
        const vector& v = values[elemI];

        if (values.size() == points.size())
        {
            writeGeometry(points, elemI, os);
        }
        else
        {
            writeGeometry(points, faces, elemI, os);
        }

        os << v[0] << ' ' << v[1] << ' ' << v[2] << nl;
    }
}


// Write sphericalTensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const sphericalTensorField& values,
    Ostream& os
) const
{
    // header
    os << "# " << fieldName;

    if (values.size() == points.size())
    {
        os  << "  POINT_DATA " << values.size()
            << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size()
            << nl;
    }

    os  << "#  ii  ";
    os << fieldName << "_ii" << endl;

    // Write data
    forAll(values, elemI)
    {
        const sphericalTensor& v = values[elemI];

        if (values.size() == points.size())
        {
            writeGeometry(points, elemI, os);
        }
        else
        {
            writeGeometry(points, faces, elemI, os);
        }

        os  << v[0] << nl;
    }
}


// Write symmTensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const symmTensorField& values,
    Ostream& os
) const
{
    // header
    os << "# " << fieldName;

    if (values.size() == points.size())
    {
        os  << "  POINT_DATA " << values.size()
            << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size()
            << nl;
    }

    os  << "#  xx  xy  xz  yy  yz ";
    for(int i=0; i<6; i++)
    {
        os << fieldName << "_" << i << "  ";
    }
    os << endl;

    // Write data
    forAll(values, elemI)
    {
        const symmTensor& v = values[elemI];

        if (values.size() == points.size())
        {
            writeGeometry(points, elemI, os);
        }
        else
        {
            writeGeometry(points, faces, elemI, os);
        }

        os  << v[0] << ' ' << v[1] << ' ' << v[2]
            << v[3] << ' ' << v[4] << ' ' << v[5]
            << nl;
    }
}


// Write tensorField in raw format
template<class Type>
void Foam::rawSurfaceWriter<Type>::writeData
(
    const fileName& fieldName,
    const pointField& points,
    const faceList& faces,
    const tensorField& values,
    Ostream& os
) const
{
    // header
    os << "# " << fieldName;

    if (values.size() == points.size())
    {
        os  << "  POINT_DATA " << values.size()
            << nl;
    }
    else
    {
        os  << "  FACE_DATA " << values.size()
            << nl;
    }

    os  << "#  xx  xy  xz  yx  yy  yz  zx  zy  zz";
    for(int i=0; i<9; i++)
    {
        os << fieldName << "_" << i << "  ";
    }
    os << endl;

    // Write data
    forAll(values, elemI)
    {
        const tensor& v = values[elemI];

        if (values.size() == points.size())
        {
            writeGeometry(points, elemI, os);
        }
        else
        {
            writeGeometry(points, faces, elemI, os);
        }

        os  << v[0] << ' ' << v[1] << ' ' << v[2]
            << v[3] << ' ' << v[4] << ' ' << v[5]
            << v[6] << ' ' << v[7] << ' ' << v[8] << nl;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSurfaceWriter<Type>::rawSurfaceWriter()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::rawSurfaceWriter<Type>::~rawSurfaceWriter()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::rawSurfaceWriter<Type>::write
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

    fileName planeFName(surfaceDir/fieldName + '_' + surfaceName + ".raw");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << planeFName << endl;
    }

    OFstream rawFile(planeFName);

    writeData(fieldName, points, faces, values, rawFile);
}


// ************************************************************************* //
