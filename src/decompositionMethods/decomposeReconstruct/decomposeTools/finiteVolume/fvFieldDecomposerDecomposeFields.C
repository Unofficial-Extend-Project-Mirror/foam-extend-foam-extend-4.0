/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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

#include "fvFieldDecomposer.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvPatchField, volMesh>& field
) const
{
    // Create and map the internal field values
    Field<Type> internalField(field.internalField(), cellAddressing_);

    // Create and map the patch field values
    PtrList<fvPatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchI]],
                    procMesh_.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null(),
                    *patchFieldDecomposerPtrs_[patchI]
                )
            );
        }
        else
        {
            // This will be a processor or passiveProcessor field,
            // depending on the patch type.  Construct it first and set the
            // data by hand
            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                // new processorFvPatchField<Type>
                (
                    procMesh_.boundary()[patchI].type(),
                    procMesh_.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );

            patchFields[patchI] ==
                Field<Type>
                (
                    field.internalField(),
                    *processorVolPatchFieldDecomposerPtrs_[patchI]
                );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, fvPatchField, volMesh> >
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
fvFieldDecomposer::decomposeField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& field
) const
{
    labelList mapAddr
    (
        labelList::subList
        (
            faceAddressing_,
            procMesh_.nInternalFaces()
        )
    );
    forAll (mapAddr, i)
    {
        mapAddr[i] -= 1;
    }

    // Create and map the internal field values
    Field<Type> internalField
    (
        field.internalField(),
        mapAddr
    );

    // Problem with addressing when a processor patch picks up both internal
    // faces and faces from cyclic boundaries. This is a bit of a hack, but
    // I cannot find a better solution without making the internal storage
    // mechanism for surfaceFields correspond to the one of faces in polyMesh
    // (i.e. using slices). HJ, date unknown, before 2004

    // Copy all face data into a single array across all patches
    Field<Type> allFaceField(field.mesh().nFaces());

    forAll (field.internalField(), i)
    {
        allFaceField[i] = field.internalField()[i];
    }

    forAll (field.boundaryField(), patchI)
    {
        const Field<Type> & p = field.boundaryField()[patchI];

        const label patchStart = field.mesh().boundaryMesh()[patchI].start();

        // Bugfix
        // Special patches soch as overset and immersed boundary have
        // their faces and values in a separate list; the list of
        // values does not correspond to faces.  Skip such patches.
        // HJ, 4/Jul/2019
        if (p.size() == field.mesh().boundaryMesh()[patchI].size())
        {
            forAll (p, i)
            {
                allFaceField[patchStart + i] = p[i];
            }
        }
    }

    // Create and map the patch field values
    PtrList<fvsPatchField<Type> > patchFields(boundaryAddressing_.size());

    forAll (boundaryAddressing_, patchI)
    {
        if (boundaryAddressing_[patchI] >= 0)
        {
            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    field.boundaryField()[boundaryAddressing_[patchI]],
                    procMesh_.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null(),
                    *patchFieldDecomposerPtrs_[patchI]
                )
            );
        }
        else
        {
            // This will be a processor or passiveProcessor field,
            // depending on the patch type.  Construct it first and set the
            // data by hand
            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                // new processorFvsPatchField<Type>
                (
                    procMesh_.boundary()[patchI].type(),
                    procMesh_.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );

            patchFields[patchI] ==
                Field<Type>
                (
                    allFaceField,
                    *processorSurfacePatchFieldDecomposerPtrs_[patchI]
                );
        }
    }

    // Create the field for the processor
    return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                field.name(),
                procMesh_.time().timeName(),
                procMesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            procMesh_,
            field.dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class GeoField>
void fvFieldDecomposer::decomposeFields
(
    const PtrList<GeoField>& fields
) const
{
    forAll (fields, fieldI)
    {
        decomposeField(fields[fieldI])().write();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
