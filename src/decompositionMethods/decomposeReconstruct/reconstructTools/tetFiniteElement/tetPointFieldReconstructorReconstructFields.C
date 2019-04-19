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

#include "tetPointFieldReconstructor.H"
#include "PtrList.H"
#include "tetPolyPatchFields.H"
#include "tetFemMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void tetPointFieldReconstructor::reconstructTetPointField
(
    GeometricField<Type, tetPolyPatchField, tetPointMesh>& reconField,
    const PtrList<GeometricField<Type, tetPolyPatchField, tetPointMesh> >&
        procFields
) const
{
    // Get references to internal and boundary field
    Field<Type>& iField = reconField.internalField();

    typename GeometricField<Type, tetPolyPatchField, tetPointMesh>::
        GeometricBoundaryField& bouField = reconField.boundaryField();

    forAll (procMeshes_, procI)
    {
        if (procFields.set(procI))
        {
            const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
                procField = procFields[procI];

            // Get processor-to-global addressing for use in rmap
            labelList procToGlobalAddr = procAddressing(procI);

            // Set the cell values in the reconstructed field
            iField.rmap
            (
                procField.internalField(),
                procToGlobalAddr
            );

            // Set the boundary patch values in the reconstructed field
            forAll(boundaryProcAddressing_[procI], patchI)
            {
                // Get patch index of the original patch
                const label curBPatch = boundaryProcAddressing_[procI][patchI];

                // check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    // If the field stores values, do the rmap
                    if (bouField[curBPatch].storesFieldData())
                    {
                        bouField[curBPatch].rmap
                        (
                            procField.boundaryField()[patchI],
                            procPatchAddressing
                            (
                                procToGlobalAddr,
                                procI,
                                patchI
                            )
                        );
                    }
                }
            }
        }
    }
}


template<class Type>
tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
tetPointFieldReconstructor::reconstructTetPointField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, tetPolyPatchField, tetPointMesh> > procFields
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
                new GeometricField<Type, tetPolyPatchField, tetPointMesh>
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
    PtrList<tetPolyPatchField<Type> > patchFields(mesh_.boundary().size());

    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            // Set the boundary patch values in the reconstructed field
            forAll(boundaryProcAddressing_[procI], patchI)
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
                            tetPolyPatchField<Type>::New
                            (
                                procFields[procI].boundaryField()[patchI],
                                mesh_.boundary()[curBPatch],
                                DimensionedField<Type, tetPointMesh>::null(),
                                tetPolyPatchFieldReconstructor
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
                tetPolyPatchField<Type>::New
                (
                    mesh_.boundary()[patchI].type(),
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, tetPointMesh>::null()
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
    tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> > treconField
    (
        new GeometricField<Type, tetPolyPatchField, tetPointMesh>
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
            Field<Type>(mesh_.nPoints()),
            patchFields
        )
    );

    // Reconstruct field
    reconstructTetPointField
    (
        treconField(),
        procFields
    );

    return treconField;
}


template<class Type>
void tetPointFieldReconstructor::reconstructElementField
(
    GeometricField<Type, elementPatchField, elementMesh>& reconField,
    const PtrList<GeometricField<Type, elementPatchField, elementMesh> >&
        procFields
) const
{
    // Get references to internal and boundary field
    Field<Type>& iField = reconField.internalField();

    typename GeometricField<Type, elementPatchField, elementMesh>::
        GeometricBoundaryField& bouField = reconField.boundaryField();

    forAll (procMeshes_, procI)
    {
        if (procFields.set(procI))
        {
            const GeometricField<Type, elementPatchField, elementMesh>&
                procField = procFields[procI];

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
                    procMeshes_[procI]().boundaryMesh()[patchI].patchSlice
                    (
                        faceProcAddressing_[procI]
                    );

                // Check if the boundary patch is not a processor patch
                if (curBPatch >= 0)
                {
                    // If the field stores values, do the rmap
                    if (bouField[curBPatch].storesFieldData())
                    {
                        const label curPatchStart =
                            mesh_().boundaryMesh()[curBPatch].start();

                        labelList reverseAddressing(cp.size());

                        forAll (cp, faceI)
                        {
                            // Subtract one to take into account offsets for
                            // face direction.
                            reverseAddressing[faceI] = cp[faceI] - 1
                                - curPatchStart;
                        }

                        bouField[curBPatch].rmap
                        (
                            procField.boundaryField()[patchI],
                            reverseAddressing
                        );
                    }
                }
            }
        }
    }
}


template<class Type>
tmp<GeometricField<Type, elementPatchField, elementMesh> >
tetPointFieldReconstructor::reconstructElementField
(
    const IOobject& fieldIoObject
) const
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, elementPatchField, elementMesh> > procFields
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
                new GeometricField<Type, elementPatchField, elementMesh>
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
    PtrList<elementPatchField<Type> > patchFields(mesh_.boundary().size());


    forAll (procMeshes_, procI)
    {
        if (procMeshes_.set(procI))
        {
            // Set the boundary patch values in the reconstructed field
            forAll(boundaryProcAddressing_[procI], patchI)
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
                            elementPatchField<Type>::New
                            (
                                procFields[procI].boundaryField()[patchI],
                                mesh_.boundary()[curBPatch],
                                DimensionedField<Type, elementMesh>::null(),
                                tetPolyPatchFieldReconstructor
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
                elementPatchField<Type>::New
                (
                    mesh_.boundary()[patchI].type(),
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, elementMesh>::null()
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
    tmp<GeometricField<Type, elementPatchField, elementMesh> > treconField
    (
        new GeometricField<Type, elementPatchField, elementMesh>
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
            Field<Type>(mesh_.nCells()),
            patchFields
        )
    );

    // Reconstruct field
    reconstructElementField
    (
        treconField(),
        procFields
    );

    return treconField;
}


template<class Type>
void tetPointFieldReconstructor::reconstructTetPointFields
(
    const IOobjectList& objects
) const
{
    word fieldClassName
    (
        GeometricField<Type, tetPolyPatchField, tetPointMesh>::typeName
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

            reconstructTetPointField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}


template<class Type>
void tetPointFieldReconstructor::reconstructElementFields
(
    const IOobjectList& objects
) const
{
    word fieldClassName
    (
        GeometricField<Type, elementPatchField, elementMesh>::typeName
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

            reconstructElementField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
