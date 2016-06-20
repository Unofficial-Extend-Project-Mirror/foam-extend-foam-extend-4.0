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

#include "solidContactFvPatchVectorField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::solidContactFvPatchVectorField::
zoneField
(
    const label zoneIndex,
    const label patchIndex,
    const Field<Type>& patchField
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faceZone& fZone = mesh.faceZones()[zoneIndex];

    tmp<Field<Type> > tZoneField
    (
        new Field<Type>(fZone.size(), pTraits<Type>::zero)
    );
    Field<Type>& zField = tZoneField();

    if (!Pstream::parRun())
    {
        zField = patchField;
    }
    else
    {
        const label patchStart
            = mesh.boundaryMesh()[patchIndex].start();

        // put local patchField into global zoneField
        forAll(patchField, i)
        {
            zField[fZone.whichFace(patchStart + i)] = patchField[i];
        }

        reduce(zField, sumOp<Field<Type> >());
    }

    return tZoneField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::solidContactFvPatchVectorField::
patchField
(
    const label patchIndex,
    const label zoneIndex,
    const Field<Type>& zoneField
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();

    const faceZone& fZone = mesh.faceZones()[zoneIndex];

    tmp<Field<Type> > tPatchField
    (
        new Field<Type>
        (
            mesh.boundaryMesh()[patchIndex].size(),
            pTraits<Type>::zero
        )
    );
    Field<Type>& pField = tPatchField();

    if (!Pstream::parRun())
    {
        pField = zoneField;
    }
    else
    {
        const label patchStart
            = mesh.boundaryMesh()[patchIndex].start();

        // take local patchField from the global zoneField
        forAll(pField, i)
        {
            pField[i] =
                zoneField[fZone.whichFace(patchStart + i)];
        }
    }

    return tPatchField;
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::solidContactFvPatchVectorField::
zoneFaceToPointInterpolate
(
    const label zoneID,
    const Field<Type>& zoneField
) const
{
    if (zoneID == zoneIndex())
    {
        PrimitivePatchInterpolation<primitiveFacePatch> zoneInterp(zone());

        tmp<Field<Type> > tZonePointField
        (
            new Field<Type>(zone().localPoints().size(), pTraits<Type>::zero)
        );
        Field<Type>& zonePointField = tZonePointField();

        zonePointField = zoneInterp.faceToPointInterpolate(zoneField);

        return tZonePointField;
    }
    else if (zoneID == shadowZoneIndex())
    {
        PrimitivePatchInterpolation<primitiveFacePatch> zoneInterp
            (
                shadowZone()
            );

        tmp<Field<Type> > tZonePointField
        (
            new Field<Type>
            (
                shadowZone().localPoints().size(), pTraits<Type>::zero
            )
        );
        Field<Type>& zonePointField = tZonePointField();

        zonePointField = zoneInterp.faceToPointInterpolate(zoneField);

        return tZonePointField;
    }
    else
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> > Foam::"
            "solidContactFvPatchVectorField::\n"
            "zoneFaceToPointInterpolate\n"
            "("
            "    const label zoneID,\n"
            "    const Field<Type>& zoneField\n"
            ") const\n"
        )   << "zone must be master or slave zone" << abort(FatalError);
    }

    return tmp<Field<Type> >();
}


template<class Type>
Foam::tmp<Foam::Field<Type> > Foam::solidContactFvPatchVectorField::
zonePointToFaceInterpolate
(
    const label zoneID,
    const Field<Type>& zonePointField
) const
{
    if (zoneID == zoneIndex())
    {
        PrimitivePatchInterpolation<primitiveFacePatch> zoneInterp(zone());

        tmp<Field<Type> > tZoneField
        (
            new Field<Type>(zone().size(), pTraits<Type>::zero)
        );
        Field<Type>& zoneField = tZoneField();

        zoneField = zoneInterp.pointToFaceInterpolate(zonePointField);

        return tZoneField;
    }
    else if (zoneID == shadowZoneIndex())
    {
        PrimitivePatchInterpolation<primitiveFacePatch> zoneInterp
            (
                shadowZone()
            );

        tmp<Field<Type> > tZoneField
        (
            new Field<Type>(shadowZone().size(), pTraits<Type>::zero)
        );
        Field<Type>& zoneField = tZoneField();

        zoneField = zoneInterp.pointToFaceInterpolate(zonePointField);

        return tZoneField;
    }
    else
    {
        FatalErrorIn
        (
            "template<class Type>\n"
            "Foam::tmp<Foam::Field<Type> > Foam::"
            "solidContactFvPatchVectorField::\n"
            "zonePointToFaceInterpolate\n"
            "("
            "    const label zoneID,\n"
            "    const Field<Type>& zonePointField\n"
            ") const\n"
        )   << "zone must be master or slave zone" << abort(FatalError);
    }

    return tmp<Field<Type> >();
}


// ************************************************************************* //
