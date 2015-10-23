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

Class
    volSurfaceMapping

Description
    Volume to surface and surface to volume mapping

Author
    Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "volSurfaceMapping.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<template<class> class PatchField, class Type>
Foam::tmp<Foam::Field<Type> > Foam::volSurfaceMapping::mapToSurface
(
    const FieldField<PatchField, Type>& df
) const
{
    // Grab labels for all faces in faMesh
    const labelList& faceLabels = mesh_.faceLabels();

    tmp<Field<Type> > tresult
    (
        new Field<Type>
        (
            faceLabels.size(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& result = tresult();

    // Get reference to volume mesh
    const polyMesh& pMesh = mesh_();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();

    label patchID, faceID;

    // Grab droplet cloud source by identifying patch and face
    forAll (faceLabels, i)
    {
        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceLabels[i] < pMesh.nFaces())
        {
            patchID = bm.whichPatch(faceLabels[i]);
            faceID = bm[patchID].whichFace(faceLabels[i]);

            result[i] = df[patchID][faceID];
        }
    }

    return tresult;
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const GeometricField<Type, faPatchField, areaMesh>& af,
    Foam::FieldField<Foam::fvPatchField, Type>& bf
) const
{
    // Grab labels for all faces in faMesh
    const labelList& faceLabels = mesh_.faceLabels();

    // Get reference to volume mesh
    const polyMesh& pMesh = mesh_();
    const polyBoundaryMesh& bm = pMesh.boundaryMesh();

    label patchID, faceID;

    const Field<Type>& afi = af.internalField();

    forAll (faceLabels, i)
    {
        // Escape if face is beyond active faces, eg belongs to a face zone
        if (faceLabels[i] < pMesh.nFaces())
        {
            patchID = bm.whichPatch(faceLabels[i]);
            faceID = bm[patchID].whichFace(faceLabels[i]);

            bf[patchID][faceID] = afi[i];
        }
    }
}


template<class Type>
void Foam::volSurfaceMapping::mapToVolume
(
    const tmp<GeometricField<Type, faPatchField, areaMesh> >& taf,
    Foam::FieldField<Foam::fvPatchField, Type>& bf
) const
{
    mapToVolume(taf(), bf);

    taf.clear();
}


// ************************************************************************* //
