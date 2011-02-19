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

#include "faFieldReconstructor.H"
#include "Time.H"
#include "PtrList.H"
#include "faPatchFields.H"
#include "emptyFaPatch.H"
#include "emptyFaPatchField.H"
#include "emptyFaePatchField.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh> >
Foam::faFieldReconstructor::reconstructFaAreaField
(
    const IOobject& fieldIoObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, faPatchField, areaMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
    {
        procFields.set
        (
            procI,
            new GeometricField<Type, faPatchField, areaMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[procI].time().timeName(),
                    procMeshes_[procI](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[procI]
            )
        );
    }

    // Create the internalField
    Field<Type> internalField(mesh_.nFaces());

    // Create the patch fields
    PtrList<faPatchField<Type> > patchFields(mesh_.boundary().size());


    // Create global mesh patchs starts

    labelList gStarts(mesh_.boundary().size(), -1);

    if (mesh_.boundary().size() > 0)
    {
        gStarts[0] = mesh_.nInternalEdges();
    }

    for(label i=1; i<mesh_.boundary().size(); i++)
    {
        gStarts[i] = gStarts[i-1] + mesh_.boundary()[i-1].labelList::size();
    }

    forAll (procMeshes_, procI)
    {
        const GeometricField<Type, faPatchField, areaMesh>& procField =
            procFields[procI];

        // Set the face values in the reconstructed field
        internalField.rmap
        (
            procField.internalField(),
            faceProcAddressing_[procI]
        );

        

        // Set the boundary patch values in the reconstructed field

        labelList starts(procMeshes_[procI].boundary().size(), -1);

        if(procMeshes_[procI].boundary().size() > 0)
        {
            starts[0] = procMeshes_[procI].nInternalEdges();
        }

        for(label i=1; i<procMeshes_[procI].boundary().size(); i++)
        {
            starts[i] = 
                starts[i-1] 
              + procMeshes_[procI].boundary()[i-1].labelList::size();
        }

        forAll(boundaryProcAddressing_[procI], patchI)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[procI][patchI];

            // Get addressing slice for this patch

//             const labelList::subList cp =
//                 procMeshes_[procI].boundary()[patchI].patchSlice
//                 (
//                     edgeProcAddressing_[procI]
//                 );

            const labelList::subList cp =
                labelList::subList
                (
                    edgeProcAddressing_[procI], 
                    procMeshes_[procI].boundary()[patchI].size(), 
                    starts[patchI]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        faPatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, areaMesh>::null(),
                            faPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                const label curPatchStart = gStarts[curBPatch];
//                     mesh_.boundary()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.  
//                     reverseAddressing[edgeI] = cp[edgeI] - 1 - curPatchStart;
                    reverseAddressing[edgeI] = cp[edgeI] - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchI],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchI];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.  
//                     label curE = cp[edgeI] - 1;
                    label curE = cp[edgeI];

                    // Is the face on the boundary?
                    if (curE >= mesh_.nInternalEdges())
                    {
//                         label curBPatch = mesh_.boundary().whichPatch(curE);
                        label curBPatch = -1;

                        forAll (mesh_.boundary(), pI)
                        {
                            if
                            (
                                curE >= gStarts[pI]
                             && curE < 
                                (
                                    gStarts[pI] 
                                  + mesh_.boundary()[pI].labelList::size()
                                )
                            )
                            {
                                curBPatch = pI;
                            }
                        }

                        if (!patchFields(curBPatch))
                        {
                            patchFields.set
                            (
                                curBPatch,
                                faPatchField<Type>::New
                                (
                                    mesh_.boundary()[curBPatch].type(),
                                    mesh_.boundary()[curBPatch],
                                    DimensionedField<Type, areaMesh>::null()
                                )
                            );
                        }

                        // add the edge
//                         label curPatchEdge =
//                             mesh_.boundary()
//                                 [curBPatch].whichEdge(curE);

                        label curPatchEdge = curE - gStarts[curBPatch];

                        patchFields[curBPatch][curPatchEdge] =
                            curProcPatch[edgeI];
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchI)
    {
        // add empty patches
        if
        (
            typeid(mesh_.boundary()[patchI]) == typeid(emptyFaPatch)
         && !patchFields(patchI)
        )
        {
            patchFields.set
            (
                patchI,
                faPatchField<Type>::New
                (
                    emptyFaPatchField<Type>::typeName,
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, areaMesh>::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, faPatchField, areaMesh> >
    (
        new GeometricField<Type, faPatchField, areaMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_.time().timeName(),
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
Foam::tmp<Foam::GeometricField<Type, Foam::faePatchField, Foam::edgeMesh> >
Foam::faFieldReconstructor::reconstructFaEdgeField
(
    const IOobject& fieldIoObject
)
{
    // Read the field for all the processors
    PtrList<GeometricField<Type, faePatchField, edgeMesh> > procFields
    (
        procMeshes_.size()
    );

    forAll (procMeshes_, procI)
    {
        procFields.set
        (
            procI,
            new GeometricField<Type, faePatchField, edgeMesh>
            (
                IOobject
                (
                    fieldIoObject.name(),
                    procMeshes_[procI].time().timeName(),
                    procMeshes_[procI](),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                ),
                procMeshes_[procI]
            )
        );
    }


    // Create the internalField
    Field<Type> internalField(mesh_.nInternalEdges());

    // Create the patch fields
    PtrList<faePatchField<Type> > patchFields(mesh_.boundary().size());


    labelList gStarts(mesh_.boundary().size(), -1);

    if(mesh_.boundary().size() > 0)
    {
        gStarts[0] = mesh_.nInternalEdges();
    }

    for(label i=1; i<mesh_.boundary().size(); i++)
    {
        gStarts[i] = gStarts[i-1] + mesh_.boundary()[i-1].labelList::size();
    }


    forAll (procMeshes_, procI)
    {
        const GeometricField<Type, faePatchField, edgeMesh>& procField =
            procFields[procI];

        // Set the face values in the reconstructed field

        // It is necessary to create a copy of the addressing array to
        // take care of the face direction offset trick.
        // 
        {
            labelList curAddr(edgeProcAddressing_[procI]);

//             forAll (curAddr, addrI)
//             {
//                 curAddr[addrI] -= 1;
//             }

            internalField.rmap
            (
                procField.internalField(),
                curAddr
            );
        }

        // Set the boundary patch values in the reconstructed field

        labelList starts(procMeshes_[procI].boundary().size(), -1);

        if(procMeshes_[procI].boundary().size() > 0)
        {
            starts[0] = procMeshes_[procI].nInternalEdges();
        }

        for(label i=1; i<procMeshes_[procI].boundary().size(); i++)
        {
            starts[i] = 
                starts[i-1] 
              + procMeshes_[procI].boundary()[i-1].labelList::size();
        }

        forAll(boundaryProcAddressing_[procI], patchI)
        {
            // Get patch index of the original patch
            const label curBPatch = boundaryProcAddressing_[procI][patchI];

            // Get addressing slice for this patch

//             const labelList::subList cp =
//                 procMeshes_[procI].boundary()[patchI].patchSlice
//                 (
//                     faceProcAddressing_[procI]
//                 );

            const labelList::subList cp =
                labelList::subList
                (
                    edgeProcAddressing_[procI], 
                    procMeshes_[procI].boundary()[patchI].size(), 
                    starts[patchI]
                );

            // check if the boundary patch is not a processor patch
            if (curBPatch >= 0)
            {
                // Regular patch. Fast looping

                if (!patchFields(curBPatch))
                {
                    patchFields.set
                    (
                        curBPatch,
                        faePatchField<Type>::New
                        (
                            procField.boundaryField()[patchI],
                            mesh_.boundary()[curBPatch],
                            DimensionedField<Type, edgeMesh>::null(),
                            faPatchFieldReconstructor
                            (
                                mesh_.boundary()[curBPatch].size(),
                                procField.boundaryField()[patchI].size()
                            )
                        )
                    );
                }

                const label curPatchStart = gStarts[curBPatch];
//                     mesh_.boundary()[curBPatch].start();

                labelList reverseAddressing(cp.size());

                forAll(cp, edgeI)
                {
                    // Subtract one to take into account offsets for
                    // face direction.  
//                     reverseAddressing[faceI] = cp[faceI] - 1 - curPatchStart;
                    reverseAddressing[edgeI] = cp[edgeI] - curPatchStart;
                }

                patchFields[curBPatch].rmap
                (
                    procField.boundaryField()[patchI],
                    reverseAddressing
                );
            }
            else
            {
                const Field<Type>& curProcPatch =
                    procField.boundaryField()[patchI];

                // In processor patches, there's a mix of internal faces (some
                // of them turned) and possible cyclics. Slow loop
                forAll(cp, edgeI)
                {
//                     label curF = cp[edgeI] - 1;
                    label curE = cp[edgeI];

                    // Is the face turned the right side round
                    if (curE >= 0)
                    {
                        // Is the face on the boundary?
                        if (curE >= mesh_.nInternalEdges())
                        {
//                             label curBPatch =
//                                 mesh_.boundary().whichPatch(curF);

                            label curBPatch = -1;

                            forAll (mesh_.boundary(), pI)
                            {
                                if
                                (
                                    curE >= gStarts[pI]
                                 && curE < 
                                    (
                                        gStarts[pI] 
                                      + mesh_.boundary()[pI].labelList::size()
                                    )
                                )
                                {
                                    curBPatch = pI;
                                }
                            }

                            if (!patchFields(curBPatch))
                            {
                                patchFields.set
                                (
                                    curBPatch,
                                    faePatchField<Type>::New
                                    (
                                        mesh_.boundary()[curBPatch].type(),
                                        mesh_.boundary()[curBPatch],
                                        DimensionedField<Type, edgeMesh>
                                           ::null()
                                    )
                                );
                            }

                            // add the face
//                             label curPatchFace =
//                                 mesh_.boundary()
//                                 [curBPatch].whichEdge(curF);

                            label curPatchEdge = curE - gStarts[curBPatch];

                            patchFields[curBPatch][curPatchEdge] =
                                curProcPatch[edgeI];
                        }
                        else
                        {
                            // Internal face
                            internalField[curE] = curProcPatch[edgeI];
                        }
                    }
                }
            }
        }
    }

    forAll(mesh_.boundary(), patchI)
    {
        // add empty patches
        if
        (
            typeid(mesh_.boundary()[patchI]) == typeid(emptyFaPatch)
         && !patchFields(patchI)
        )
        {
            patchFields.set
            (
                patchI,
                faePatchField<Type>::New
                (
                    emptyFaePatchField<Type>::typeName,
                    mesh_.boundary()[patchI],
                    DimensionedField<Type, edgeMesh>::null()
                )
            );
        }
    }


    // Now construct and write the field
    // setting the internalField and patchFields
    return tmp<GeometricField<Type, faePatchField, edgeMesh> >
    (
        new GeometricField<Type, faePatchField, edgeMesh>
        (
            IOobject
            (
                fieldIoObject.name(),
                mesh_.time().timeName(),
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


// Reconstruct and write all area fields
template<class Type>
void Foam::faFieldReconstructor::reconstructFaAreaFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName =
        GeometricField<Type, faPatchField, areaMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        for
        (
            IOobjectList::const_iterator fieldIter = fields.begin();
            fieldIter != fields.end();
            ++fieldIter
        )
        {
            Info << "        " << fieldIter()->name() << endl;

            reconstructFaAreaField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}

// Reconstruct and write all edge fields
template<class Type>
void Foam::faFieldReconstructor::reconstructFaEdgeFields
(
    const IOobjectList& objects
)
{
    const word& fieldClassName =
        GeometricField<Type, faePatchField, edgeMesh>::typeName;

    IOobjectList fields = objects.lookupClass(fieldClassName);

    if (fields.size())
    {
        Info<< "    Reconstructing " << fieldClassName << "s\n" << endl;

        for
        (
            IOobjectList::const_iterator fieldIter = fields.begin();
            fieldIter != fields.end();
            ++fieldIter
        )
        {
            Info<< "        " << fieldIter()->name() << endl;

            reconstructFaEdgeField<Type>(*fieldIter())().write();
        }

        Info<< endl;
    }
}


// ************************************************************************* //
