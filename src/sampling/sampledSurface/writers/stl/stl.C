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

#include "stl.H"
#include "fileName.H"
#include "OFstream.H"
#include "faceList.H"
#include "OSspecific.H"
#include "triSurface.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::stl<Type>::stl()
:
    surfaceWriter<Type>()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::stl<Type>::~stl()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::stl<Type>::write
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

    fileName planeFName(surfaceDir/fieldName + '_' + surfaceName + ".stl");

    if (verbose)
    {
        Info<< "Writing field " << fieldName << " to " << planeFName << endl;
    }

    // Convert faces to triangles.
    DynamicList<labelledTri> tris(faces.size());

    forAll(faces, i)
    {
        const face& f = faces[i];

        faceList triFaces(f.nTriangles(points));
        label nTris = 0;
        f.triangles(points, nTris, triFaces);

        forAll(triFaces, triI)
        {
            const face& tri = triFaces[triI];
            tris.append(labelledTri(tri[0], tri[1], tri[2], 0));
        }
    }

    triSurface
    (
        tris.shrink(),
        geometricSurfacePatchList
        (
            1,
            geometricSurfacePatch
            (
                "patch",                            // geometricType
                string::validate<word>(fieldName),  // fieldName
                0                                   // index
            )
        ),
        points
    ).write(planeFName);
}


// ************************************************************************* //
