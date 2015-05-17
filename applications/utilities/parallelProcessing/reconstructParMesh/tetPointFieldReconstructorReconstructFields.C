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
tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
tetPointFieldReconstructor::reconstructTetPointField
(
    const IOobject& fieldIoObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, tetPolyPatchField, tetPointMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
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


    // Create the internalField
    Field<Type> internalField(mesh_.nPoints());

    // Create the patch fields
    PtrList<tetPolyPatchField<Type> > patchFields(mesh_.boundary().size());


    forAll (procMeshes_, procI)
    {
        const GeometricField<Type, tetPolyPatchField, tetPointMesh>&
            procField = procFields[procI];

        // Get processor-to-global addressing for use in rmap
        labelList procToGlobalAddr = procAddressing(procI);

        // Set the cell values in the reconstructed field
        internalField.rmap
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
                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        tetPolyPatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, tetPointMesh>::null(),
                            tetPolyPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                // If the field stores values, do the rmap
                if (patchFields[curBPatch].storesFieldData())
                {
                    patchFields[curBPatch].rmap
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

    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, tetPolyPatchField, tetPointMesh> >
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
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
tmp<GeometricField<Type, elementPatchField, elementMesh> >
tetPointFieldReconstructor::reconstructElementField
(
    const IOobject& fieldIoObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, elementPatchField, elementMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
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


    // Create the internalField
    Field<Type> internalField(mesh_.nCells());

    // Create the patch fields
    PtrList<elementPatchField<Type> > patchFields(mesh_.boundary().size());


    forAll (procMeshes_, procI)
    {
        const GeometricField<Type, elementPatchField, elementMesh>&
            procField = procFields[procI];

        // Set the cell values in the reconstructed field
        internalField.rmap
        (
            procField.internalField(),
            cellProcAddressing_[procI]
        );

        // Set the boundary patch values in the reconstructed field
        forAll(boundaryProcAddressing_[procI], patchI)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[procI][patchI];

            // Get addressing slice for this patch
            const labelList::subList cp =
                procMeshes_[procI]().boundaryMesh()[patchI].patchSlice
                (
                    faceProcAddressing_[procI]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        elementPatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, elementMesh>::null(),
                            tetPolyPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                // If the field stores values, do the rmap
                if (patchFields[curBPatch].storesFieldData())
                {
                    const label curPatchStart =
                        mesh_().boundaryMesh()[curBPatch].start();

                    labelList reverseAddressing(cp.size());

                    forAll(cp, faceI)
                    {
                        // Subtract one to take into account offsets for
                        // face direction.
                        reverseAddressing[faceI] = cp[faceI] - 1
                            - curPatchStart;
                    }

                    patchFields[curBPatch].rmap
                    (
                        procField.boundaryField()[patchI],
                        reverseAddressing
                    );
                }
            }
        }
    }

    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, elementPatchField, elementMesh> >
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
            procFields[0].dimensions(),
            internalField,
            patchFields
        )
    );
}


template<class Type>
void tetPointFieldReconstructor::reconstructTetPointFields
(
    const IOobjectList& objects
)
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
)
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
