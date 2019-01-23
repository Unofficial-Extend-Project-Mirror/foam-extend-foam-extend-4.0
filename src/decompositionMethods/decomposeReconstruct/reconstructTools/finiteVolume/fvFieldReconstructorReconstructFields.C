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

#include "fvFieldReconstructor.H"
#include "foamTime.H"
#include "PtrList.H"
#include "fvPatchFields.H"
#include "emptyFvPatch.H"
#include "emptyFvPatchField.H"
#include "emptyFvsPatchField.H"
#include "processorFvPatch.H"
#include "processorFvPatchField.H"
#include "processorFvsPatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::fvFieldReconstructor::reconstructField
(
    GeometricField<Type, fvPatchField, volMesh>& reconField,
    const PtrList<GeometricField<Type, fvPatchField, volMesh> >& procFields
) const
{
    // Get references to internal and boundary field
    Field<Type>& iField = reconField;

    typename GeometricField<Type, fvPatchField, volMesh>::
        GeometricBoundaryField& bouField = reconField.boundaryFieldNoStoreOldTimes();

    forAll (procFields, procI)
    {
        if (procFields.set(procI))
        {
            const GeometricField<Type, fvPatchField, volMesh>& procField =
                procFields[procI];

            // Set the cell values in the reconstructed field
            iField.rmap
            (
                procField.internalField(),
                cellProcAddressing_[procI]
            );

            // Set the boundary patch values in the reconstructed field
            forAll (boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // Get addressing slice for this patch
                const labelList::subList cp =
                    procMeshes_[procI].boundary()[patchI].patchSlice
                    (
                        faceProcAddressing_[procI]
                    );

                // Check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    // Regular patch. Fast looping
                    const label curPatchStart =
                        mesh_.boundaryMesh()[curBPatch].start();

                    labelList reverseAddressing(cp.size());

                    forAll (cp, faceI)
                    {
                        // Subtract one to take into account offsets for
                        // face direction.
                        reverseAddressing[faceI] =
                            cp[faceI] - 1 - curPatchStart;
                    }

                    bouField[curBPatch].rmap
                    (
                        procField.boundaryField()[patchI],
                        reverseAddressing
                    );
                }
                else
                {
                    const Field<Type>& curProcPatch =
                        procField.boundaryField()[patchI];

                    // General mapping: patch to patch or patch to/from internal

                    // In processor patches, there's a mix of internal faces
                    // (some of them turned) and possible cyclics. Slow loop
                    forAll (cp, faceI)
                    {
                        // Subtract one to take into account offsets for
                        // face direction.
                        label curF = cp[faceI] - 1;

                        // Is the face on the boundary?
                        if (curF >= mesh_.nInternalFaces())
                        {
                            label curBPatch =
                                mesh_.boundaryMesh().whichPatch(curF);

                            // Add the face
                            label curPatchFace =
                                mesh_.boundaryMesh()[curBPatch].whichFace(curF);

                            bouField[curBPatch][curPatchFace] =
                                curProcPatch[faceI];
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
void Foam::fvFieldReconstructor::reconstructField
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& reconField,
    const PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> >& procFields
) const
{
    // Create the internalField
    Field<Type>& iField = reconField;

    typename GeometricField<Type, fvsPatchField, surfaceMesh>::
        GeometricBoundaryField& bouField = reconField.boundaryFieldNoStoreOldTimes();

    forAll (procFields, procI)
    {
        if (procFields.set(procI))
        {
            const GeometricField<Type, fvsPatchField, surfaceMesh>& procField =
                procFields[procI];

            // Set the face values in the reconstructed field

            // It is necessary to create a copy of the addressing array to
            // take care of the face direction offset trick.
            //
            {
                labelList curAddr(faceProcAddressing_[procI]);

                forAll (curAddr, addrI)
                {
                    curAddr[addrI] -= 1;
                }

                iField.rmap
                (
                    procField.internalField(),
                    curAddr
                );
            }

            // Set the boundary patch values in the reconstructed field
            forAll (boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // Get addressing slice for this patch
                const labelList::subList cp =
                    procMeshes_[procI].boundary()[patchI].patchSlice
                    (
                        faceProcAddressing_[procI]
                    );

                // Check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    // Regular patch. Fast looping
                    const label curPatchStart =
                        mesh_.boundaryMesh()[curBPatch].start();

                    labelList reverseAddressing(cp.size());

                    forAll (cp, faceI)
                    {
                        // Subtract one to take into account offsets for
                        // face direction.
                        reverseAddressing[faceI] =
                            cp[faceI] - 1 - curPatchStart;
                    }

                    bouField[curBPatch].rmap
                    (
                        procField.boundaryField()[patchI],
                        reverseAddressing
                    );
                }
                else
                {
                    const Field<Type>& curProcPatch =
                        procField.boundaryField()[patchI];

                    // In processor patches, there's a mix of internal faces
                    // (some of them turned) and possible cyclics. Slow loop
                    forAll (cp, faceI)
                    {
                        label curF = cp[faceI] - 1;

                        // Is the face turned the right side round
                        if (curF >= 0)
                        {
                            // Is the face on the boundary?
                            if (curF >= mesh_.nInternalFaces())
                            {
                                label curBPatch =
                                    mesh_.boundaryMesh().whichPatch(curF);

                                // Add the face
                                label curPatchFace =
                                    mesh_.boundaryMesh()
                                    [curBPatch].whichFace(curF);

                                bouField[curBPatch][curPatchFace] =
                                    curProcPatch[faceI];
                            }
                            else
                            {
                                // Internal face
                                iField[curF] = curProcPatch[faceI];
                            }
                        }
                    }
                }
            }
        }
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvPatchField, Foam::volMesh> >
Foam::fvFieldReconstructor::reconstructFvVolumeField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvPatchField, volMesh> > procFields
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
                new GeometricField<Type, fvPatchField, volMesh>
                (
                    IOobject
                    (
                        fieldIoObject.name(),
                        procMeshes_[procI].time().timeName(),
                        procMeshes_[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    procMeshes_[procI]
                )
            );
        }
    }

    // Create the patch fields
    PtrList<fvPatchField<Type> > patchFields(mesh_.boundary().size());

    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            // Set the boundary patch in the reconstructed field
            forAll (boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // Check if the boundary field is a regular patch field
                if (curBPatch >= 0)
                {
                    // Regular patch. Clone with map
                    if (!patchFields(curBPatch))
                    {
                        patchFields.set
                        (
                            curBPatch,
                            fvPatchField<Type>::New
                            (
                                procFields[procI].boundaryField()[patchI],
                                mesh_.boundary()[curBPatch],
                                DimensionedField<Type, volMesh>::null(),
                                fvPatchFieldReconstructor
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
                fvPatchField<Type>::New
                (
                    mesh_.boundary()[patchI].type(),
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
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
    tmp<GeometricField<Type, fvPatchField, volMesh> > treconField
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            procFields[firstValidMesh].dimensions(),
            Field<Type>(mesh_.nCells()),
            patchFields
        )
    );

    // Reconstruct field
    this->reconstructField
    (
        treconField(),
        procFields
    );

    return treconField;
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::fvFieldReconstructor::reconstructFvSurfaceField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, fvsPatchField, surfaceMesh> > procFields
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
                new GeometricField<Type, fvsPatchField, surfaceMesh>
                (
                    IOobject
                    (
                        fieldIoObject.name(),
                        procMeshes_[procI].time().timeName(),
                        procMeshes_[procI],
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE
                    ),
                    procMeshes_[procI]
                )
            );
        }
    }

    // Create the patch fields
    PtrList<fvsPatchField<Type> > patchFields(mesh_.boundary().size());

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
                    // Regular patch. Fast looping
                    if (!patchFields(curBPatch))
                    {
                        patchFields.set
                        (
                            curBPatch,
                            fvsPatchField<Type>::New
                            (
                                procFields[procI].boundaryField()[patchI],
                                mesh_.boundary()[curBPatch],
                                DimensionedField<Type, surfaceMesh>::null(),
                                fvPatchFieldReconstructor
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
                fvsPatchField<Type>::New
                (
                    mesh_.boundary()[patchI].type(),
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
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

    // Now construct the field with patches
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > treconField
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            procFields[firstValidMesh].dimensions(),
            Field<Type>(mesh_.nInternalFaces()),
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


// Reconstruct and write all/selected volume fields
template<class Type>
void Foam::fvFieldReconstructor::reconstructFvVolumeFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        GeometricField<Type, fvPatchField, volMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                !selectedFields.size()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvVolumeField<Type>(*fieldIter())().write();
            }
        }
        Info<< endl;
    }
}


// Reconstruct and write all/selected surface fields
template<class Type>
void Foam::fvFieldReconstructor::reconstructFvSurfaceFields
(
    const IOobjectList& objects,
    const HashSet<word>& selectedFields
) const
{
    const word& fieldClassName =
        GeometricField<Type, fvsPatchField, surfaceMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        forAllConstIter(IOobjectList, fields, fieldIter)
        {
            if
            (
                !selectedFields.size()
             || selectedFields.found(fieldIter()->name())
            )
            {
                Info<< "        " << fieldIter()->name() << endl;

                reconstructFvSurfaceField<Type>(*fieldIter())().write();
            }
        }

        Info<< endl;
    }
}


// ************************************************************************* //
