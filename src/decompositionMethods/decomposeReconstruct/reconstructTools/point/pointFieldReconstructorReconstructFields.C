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

#include "pointFieldReconstructor.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::pointFieldReconstructor::reconstructField
(
    GeometricField<Type, pointPatchField, pointMesh>& reconField,
    const PtrList<GeometricField<Type, pointPatchField, pointMesh> >& procFields
) const
{
    // Get references to internal and boundary field
    Field<Type>& iField = reconField.internalField();

    typename GeometricField<Type, pointPatchField, pointMesh>::
        GeometricBoundaryField& bouField = reconField.boundaryField();

    forAll (procMeshes_, procI)
    {
        if (procFields.set(procI))
        {
            const GeometricField<Type, pointPatchField, pointMesh>&
                procField = procFields[procI];

            // Get processor-to-global addressing for use in rmap
            const labelList& procToGlobalAddr = pointProcAddressing_[procI];

            // Set the cell values in the reconstructed field
            iField.rmap
            (
                procField.internalField(),
                procToGlobalAddr
            );

            // Set the boundary patch values in the reconstructed field
            forAll (boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // Check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    bouField[curBPatch].rmap
                    (
                        procFields[procI].boundaryField()[patchI],
                        patchPointAddressing_[procI][patchI]
                    );
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::pointPatchField, Foam::pointMesh> >
Foam::pointFieldReconstructor::reconstructField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, pointPatchField, pointMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            procFields.set
            (
                procI,
                new GeometricField<Type, pointPatchField, pointMesh>
                (
                    IOobject
                    (
                        fieldIoObject.name(),
                        procMeshes_[procI]().time().timeName(),
                        procMeshes_[procI](),
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    procMeshes_[procI]
                )
            );
        }
    }

    // Create the patch fields
    PtrList<pointPatchField<Type> > patchFields(mesh_.boundary().size());

    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            // Set the boundary patch in the reconstructed field
            forAll (boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // Check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    if (!patchFields(curBPatch))
                    {
                        patchFields.set
                        (
                            curBPatch,
                            pointPatchField<Type>::New
                            (
                                procFields[procI].boundaryField()[patchI],
                                mesh_.boundary()[curBPatch],
                                DimensionedField<Type, pointMesh>::null(),
                                pointPatchFieldReconstructor
                                (
                                    mesh_.boundary()[curBPatch].size(),
                                    procFields[procI].boundaryField()
                                        [patchI].size()
                                )
                            )
                        );
                    }
                }
            }
        }
    }

    // Add missing patch fields
    forAll (mesh_.boundary(), patchI)
    {
        if (!patchFields(patchI))
        {
            patchFields.set
            (
                patchI,
                pointPatchField<Type>::New
                (
                    mesh_.boundary()[patchI].type(),
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, pointMesh>::null()
                )
            );
        }
    }

    label firstValidMesh = 0;

    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            firstValidMesh = procI;
            break;
        }
    }

    // Construct the reconstructed field with patches
    tmp<GeometricField<Type, pointPatchField, pointMesh> > treconField
    (
        new GeometricField<Type, pointPatchField, pointMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_().time().timeName(),
                mesh_(),
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            procFields[firstValidMesh].dimensions(),
            Field<Type>(mesh_.size()),
            patchFields
        )
    );

    // Reconstruct field
    reconstructField
    (
        treconField(),
        procFields
    );

    return treconField;
}


// Reconstruct and write all point fields
template<class Type>
void Foam::pointFieldReconstructor::reconstructFields
(
    const IOobjectList& objects
) const
{
    word fieldClassName
    (
        GeometricField<Type, pointPatchField, pointMesh>::typeName
    );

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        for
        (
            IOobjectList::iterator fieldIter = fields.begin();
            fieldIter != fields.end();
            ++fieldIter
        )
        {
            Info<< "        " << fieldIter()->name() << endl;

            reconstructField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}


// ************************************************************************* //
