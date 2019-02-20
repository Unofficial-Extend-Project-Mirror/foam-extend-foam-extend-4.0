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

#include "topoChangerFvMesh.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class GeoField, class Decomposer>
void Foam::topoChangerFvMesh::sendFields
(
    const HashTable<const GeoField*>& geoFields,
    const Decomposer& decomposer,
    Ostream& toProc
) const
{
    // Send number of fields
    toProc<< geoFields.size() << nl << token::BEGIN_LIST << nl;

    label fI = 0;

    forAllConstIter
    (
        typename HashTable<const GeoField*>,
        geoFields,
        iter
    )
    {
        // Encapsulate a field within a dictionary in order
        // to read it correctly on the receiving side
        toProc << iter()->name() << nl
            << token::BEGIN_BLOCK << nl
            << decomposer.decomposeField(*(iter()))
            << nl << token::END_BLOCK
            << nl;

        fI++;
    }
    toProc<< token::END_LIST << nl << nl;
}


template<class GeoField, class Decomposer>
void Foam::topoChangerFvMesh::insertFields
(
    const HashTable<const GeoField*>& geoFields,
    const Decomposer& decomposer,
    List<PtrList<GeoField> >& localFields
) const
{
    // Note:
    // In the insertion, the first index in the receivedFields is the field
    // index and the second index is the processor index.
    // In this way, receivedFields[fI] gives the PtrList of processor fields
    // to be reconstructed into a single field
    // HJ, 30/Apr/2018

    // Insertion of fields happens processor-per-processor, meaning that the
    // first index is varied and procIndex remains constant
    localFields.setSize(geoFields.size());

    label fI = 0;

    forAllConstIter
    (
        typename HashTable<const GeoField*>,
        geoFields,
        iter
    )
    {
        localFields[fI].set
        (
            Pstream::myProcNo(),
            decomposer.decomposeField(*(iter()))
        );

        fI++;
    }
}


template<class GeoMesh, class GeoField>
void Foam::topoChangerFvMesh::receiveFields
(
    const label procIndex,
    List<PtrList<GeoField> >& receivedFields,
    const GeoMesh& procMesh,
    Istream& fromProc
) const
{
    // Note:
    // In the insertion, the first index in the receivedFields is the field
    // index and the second index is the processor index.
    // In this way, receivedFields[fI] gives the PtrList of processor fields
    // to be reconstructed into a single field
    // HJ, 30/Apr/2018

    // Insertion of fields happens processor-per-processor, meaning that the
    // first index is varied and procIndex remains constant

    // Receive number of fields
    label nScalarFields = readLabel(fromProc);

    receivedFields.setSize(nScalarFields);

    fromProc.readBegin("topoChangerFvMeshReceiveFields");

    forAll (receivedFields, fI)
    {
        word fieldName(fromProc);

        dictionary dict(fromProc);

        receivedFields[fI].set
        (
            procIndex,
            new GeoField
            (
                IOobject
                (
                    fieldName,
                    procMesh.time().timeName(),
                    procMesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                procMesh,
                dict
            )
        );
    }

    fromProc.readEnd("topoChangerFvMeshReceiveFields");
}


template<class GeoField, class Reconstructor>
void Foam::topoChangerFvMesh::rebuildFields
(
    const HashTable<const GeoField*>& geoFields,
    const Reconstructor& reconstructor,
    const List<PtrList<GeoField> >& receivedFields,
    const mapPolyMesh& meshMap
) const
{
    label fieldI = 0;

    // Make an fvMesh mapper
    const fvMeshMapper mapper(*this, meshMap);

    forAllConstIter
    (
        typename HashTable<const GeoField*>,
        geoFields,
        iter
    )
    {
        GeoField& masterField = const_cast<GeoField&>(*iter());

        const PtrList<GeoField>& partFields = receivedFields[fieldI];

        if (debug)
        {
            Pout<< "Rebuilding field " << masterField.name() << endl;
        }

        // Check name match.  Note: there may be holes
        word partName;

        bool nameSet = false;
        bool nameOk = true;

        forAll (partFields, pfI)
        {
            if (partFields.set(pfI))
            {
                if (!nameSet)
                {
                    partName = partFields[pfI].name();
                    nameSet = true;
                }
                else
                {
                    if (partName != partFields[pfI].name())
                    {
                        Pout<< "BAD NAME, 1: " << partName
                            << " and " << partFields[pfI].name()
                            << endl;
                        nameOk = false;
                    }
                }
            }
        }

        if (masterField.name() != partName)
        {
            // Check that all fields have the same name
            Pout<< "BAD NAME, 2 " << partName
                << " and " << masterField.name()
                << endl;
            nameOk = false;
        }

        if (!nameOk)
        {
            FatalErrorIn
            (
                "topoChangerFvMesh::rebuildFields\n"
                "(\n"
                "    const HashTable<const GeoField*>& geoFields,\n"
                "    const Reconstructor& reconstructor,\n"
                "    const List<PtrList<GeoField> >& receivedFields,\n"
                "    const mapPolyMesh& meshMap,\n"
                ") const"
            )   << "Name mismatch when rebuilding field " << masterField.name()
                << abort(FatalError);
        }

        if
        (
            meshMap.resetPatchFlag().size()
         != masterField.mesh().boundary().size()
        )
        {
            FatalErrorIn
            (
                "topoChangerFvMesh::rebuildFields\n"
                "(\n"
                "    const HashTable<const GeoField*>& geoFields,\n"
                "    const Reconstructor& reconstructor,\n"
                "    const List<PtrList<GeoField> >& receivedFields,\n"
                "    const mapPolyMesh& meshMap,\n"
                ") const"
            )   << "Bad size of patches to replace.  Boundary: "
                << masterField.mesh().boundary().size()
                << " resetPatchFlag: " << meshMap.resetPatchFlag().size()
                << abort(FatalError);
        }

        typename GeoField::GeoMeshType fieldMesh(*this);

        // Resize the master field

        // Resize internal field
        if
        (
            masterField.size()
         != GeoField::GeoMeshType::size(masterField.mesh())
        )
        {
            if (debug)
            {
                Pout<< "Resizing internal field: old size = "
                    << masterField.size()
                    << " new size = "
                    << GeoField::GeoMeshType::size(masterField.mesh())
                    << endl;
            }

            masterField.setSize
            (
                GeoField::GeoMeshType::size(masterField.mesh())
            );
        }

        // Resize boundary (number of patches)
        typename GeoField::GeometricBoundaryField& patchFields =
            masterField.boundaryFieldNoStoreOldTimes();

        if (patchFields.size() != masterField.mesh().boundary().size())
        {
            if (debug)
            {
                Pout<< "Resizing boundary field: "
                    << masterField.mesh().boundary().size()
                    << endl;
            }

            patchFields.setSize(masterField.mesh().boundary().size());
        }

        // Resize patch fields
        forAll (patchFields, patchI)
        {
            if (meshMap.resetPatchFlag()[patchI])
            {
                // Create a new constrained patch field
                if (debug)
                {
                    Pout<< "Inserting constrained patch field for patch "
                        << masterField.mesh().boundary()[patchI].name()
                        << endl;
                }

                patchFields.set
                (
                    patchI,
                    GeoField::PatchFieldType::New
                    (
                        masterField.mesh().boundary()[patchI].type(),
                        masterField.mesh().boundary()[patchI],
                        masterField
                    )
                );

                // Set field to zero to avoid uninitialised memory
                patchFields[patchI] =
                    pTraits
                    <
                        typename GeoField::PrimitiveType
                    >::zero;
            }
            else if
            (
                patchFields[patchI].size()
             != masterField.mesh().boundary()[patchI].size()
            )
            {
                // Resize patch field
                if (debug)
                {
                    Pout<< "Resizing patch field for patch "
                        << masterField.mesh().boundary()[patchI].name()
                        << " old size: " << patchFields[patchI].size()
                        << " new size: "
                        << masterField.mesh().boundary()[patchI].size()
                        << endl;
                }

                // Reset patch field size
                patchFields[patchI].autoMap
                (
                    mapper.boundaryMap()[patchI]
                );
            }
        }

        // Rebuild the master field from components
        reconstructor.reconstructField
        (
            masterField,
            partFields
        );

        // Increment field counter
        fieldI++;

        if (debug)
        {
            Pout<< "... done" << endl;
        }
    }

    // HR 14.12.18: We create new processor boundary faces from internal
    // faces. The values on these faces could be initialised by interpolation.
    // Instead we choose to fix the values by evaluating the boundaries.
    // I tried to execute evaluateCoupled() at the end of
    // fvFieldReconstructor::reconstructField, but this fails in a strange way.
    forAllConstIter
    (
        typename HashTable<const GeoField*>,
        geoFields,
        iter
    )
    {
        GeoField& masterField = const_cast<GeoField&>(*iter());
        masterField.boundaryField().evaluateCoupled();
    }
}


// ************************************************************************* //
