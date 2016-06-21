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

#include "objectRegistry.H"
#include "foamTime.H"
#include "coupledInfo.H"
#include "emptyPolyPatch.H"
#include "processorPolyPatch.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct given mesh, coupleMap and master / slave indices
template <class MeshType>
coupledInfo<MeshType>::coupledInfo
(
    const MeshType& mesh,
    const coupleMap& cMap,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_(cMap),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


// Construct from components
template <class MeshType>
coupledInfo<MeshType>::coupledInfo
(
    const MeshType& mesh,
    const bool isTwoDMesh,
    const bool isLocal,
    const bool isSend,
    const label patchIndex,
    const label mPatch,
    const label sPatch,
    const label mfzIndex,
    const label sfzIndex
)
:
    mesh_(mesh),
    builtMaps_(false),
    map_
    (
        IOobject
        (
            "coupleMap_"
          + Foam::name(mPatch)
          + "_To_"
          + Foam::name(sPatch)
          + word(isLocal ? "_Local" : "_Proc")
          + word(isSend ? "_Send" : "_Recv"),
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            true
        ),
        isTwoDMesh,
        isLocal,
        isSend,
        patchIndex,
        mPatch,
        sPatch
    ),
    masterFaceZone_(mfzIndex),
    slaveFaceZone_(sfzIndex)
{}


//- Construct given addressing
template <class MeshType>
coupledInfo<MeshType>::subMeshMapper::subMeshMapper
(
    const coupledInfo& cInfo,
    const label patchI
)
:
    sizeBeforeMapping_(cInfo.baseMesh().boundary()[patchI].size()),
    directAddressing_
    (
        SubList<label>
        (
            cInfo.map().faceMap(),
            cInfo.subMesh().boundary()[patchI].size(),
            cInfo.subMesh().boundary()[patchI].patch().start()
        )
    )
{
    // Offset indices
    label pStart = cInfo.baseMesh().boundary()[patchI].patch().start();

    forAll(directAddressing_, faceI)
    {
        directAddressing_[faceI] -= pStart;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return a const reference to the parent mesh
template <class MeshType>
const MeshType&
coupledInfo<MeshType>::baseMesh() const
{
    return mesh_;
}


// Set a new subMesh
template <class MeshType>
void coupledInfo<MeshType>::setMesh
(
    label index,
    MeshType* mesh
)
{
    subMesh_.set(mesh);
}


// Return a reference to the subMesh
template <class MeshType>
MeshType& coupledInfo<MeshType>::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("MeshType& coupledInfo::subMesh()")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


// Return a const reference to the subMesh
template <class MeshType>
const MeshType& coupledInfo<MeshType>::subMesh() const
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("const MeshType& coupledInfo::subMesh() const")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


// Return if maps have been built
template <class MeshType>
bool coupledInfo<MeshType>::builtMaps() const
{
    return builtMaps_;
}


// Set internal state of maps as built
template <class MeshType>
void coupledInfo<MeshType>::setBuiltMaps()
{
    builtMaps_ = true;
}


// Return a reference to the coupleMap
template <class MeshType>
coupleMap& coupledInfo<MeshType>::map()
{
    return map_;
}


// Return a const reference to the coupleMap
template <class MeshType>
const coupleMap& coupledInfo<MeshType>::map() const
{
    return map_;
}


// Return the master face zone ID
template <class MeshType>
label coupledInfo<MeshType>::masterFaceZone() const
{
    return masterFaceZone_;
}


// Return the slave face zone ID
template <class MeshType>
label coupledInfo<MeshType>::slaveFaceZone() const
{
    return slaveFaceZone_;
}


// Subset geometric field
template <class MeshType>
template <class GeomField, class ZeroType>
tmp<GeomField>
coupledInfo<MeshType>::subSetField
(
    const GeomField& f,
    const ZeroType& zeroValue,
    const labelList& internalMapper
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Create and map the internal-field values
    InternalField internalField(f.internalField(), internalMapper);

    // Create and map the patch field values
    label nPatches = subMesh().boundary().size();
    PtrList<PatchFieldType> patchFields(nPatches);

    // Define patch type names, assumed to be
    // common for volume and surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    // Create dummy types for initial field creation
    forAll(patchFields, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            patchFields.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    subMesh().boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                PatchFieldType::New
                (
                    PatchFieldType::calculatedType(),
                    subMesh().boundary()[patchI],
                    DimInternalField::null()
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeomField> subFld
    (
        new GeomField
        (
            IOobject
            (
                "subField_" + f.name(),
                subMesh().time().timeName(),
                subMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh(),
            f.dimensions(),
            internalField,
            patchFields
        )
    );

    // Set correct references for patch internal fields,
    // and map values from the supplied geometric field
    GeomBdyFieldType& bf = subFld().boundaryField();

    forAll(bf, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    emptyType,
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );
        }
        else
        if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    processorType,
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField()
                )
            );

            // Avoid dealing with uninitialised values
            // by artificially assigning to zero
            bf[patchI] == zeroValue;
        }
        else
        {
            bf.set
            (
                patchI,
                PatchFieldType::New
                (
                    f.boundaryField()[patchI],
                    subMesh().boundary()[patchI],
                    subFld().dimensionedInternalField(),
                    subMeshMapper(*this, patchI)
                )
            );
        }
    }

    return subFld;
}


// Subset geometric fields from registry to output stream
template <class MeshType>
template <class GeomField, class ZeroType>
void coupledInfo<MeshType>::send
(
    const wordList& fieldNames,
    const word& fieldType,
    const ZeroType& zeroValue,
    const labelList& internalMapper,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        // Fetch object from registry
        const objectRegistry& db = mesh_.thisDb();

        const GeomField& fld = db.lookupObject<GeomField>(fieldNames[i]);

        // Subset the field
        tmp<GeomField> tsubFld = subSetField(fld, zeroValue, internalMapper);

        // Send field subset through stream
        strStream
            << fieldNames[i]
            << token::NL << token::BEGIN_BLOCK
            << tsubFld
            << token::NL << token::END_BLOCK
            << token::NL;
    }

    strStream
        << token::END_BLOCK << token::NL;
}


// Set geometric field pointers from input dictionary
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::setField
(
    const wordList& fieldNames,
    const dictionary& fieldDicts,
    const label internalSize,
    PtrList<GeomField>& fields
) const
{
    typedef typename GeomField::InternalField InternalField;
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBdyFieldType;
    typedef typename GeomField::DimensionedInternalField DimInternalField;

    // Size up the pointer list
    fields.setSize(fieldNames.size());

    // Define patch type names, assumed to be
    // common for volume and surface fields
    word emptyType(emptyPolyPatch::typeName);
    word processorType(processorPolyPatch::typeName);

    forAll(fieldNames, i)
    {
        // Create and map the patch field values
        label nPatches = subMesh().boundary().size();

        // Create field parts
        PtrList<PatchFieldType> patchFields(nPatches);

        // Read dimensions
        dimensionSet dimSet
        (
            fieldDicts.subDict(fieldNames[i]).lookup("dimensions")
        );

        // Read the internal field
        InternalField internalField
        (
            "internalField",
            fieldDicts.subDict(fieldNames[i]),
            internalSize
        );

        // Create dummy types for initial field creation
        forAll(patchFields, patchI)
        {
            if (patchI == (nPatches - 1))
            {
                // Artificially set last patch
                patchFields.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        emptyType,
                        subMesh().boundary()[patchI],
                        DimInternalField::null()
                    )
                );
            }
            else
            {
                patchFields.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        PatchFieldType::calculatedType(),
                        subMesh().boundary()[patchI],
                        DimInternalField::null()
                    )
                );
            }
        }

        // Create field with dummy patches
        fields.set
        (
            i,
            new GeomField
            (
                IOobject
                (
                    fieldNames[i],
                    subMesh().time().timeName(),
                    subMesh(),
                    IOobject::NO_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                subMesh(),
                dimSet,
                internalField,
                patchFields
            )
        );

        // Set correct references for patch internal fields,
        // and fetch values from the supplied geometric field dictionaries
        GeomBdyFieldType& bf = fields[i].boundaryField();

        forAll(bf, patchI)
        {
            if (patchI == (nPatches - 1))
            {
                // Artificially set last patch
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        emptyType,
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            if (isA<processorPolyPatch>(subMesh().boundary()[patchI].patch()))
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        processorType,
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField()
                    )
                );
            }
            else
            {
                bf.set
                (
                    patchI,
                    PatchFieldType::New
                    (
                        subMesh().boundary()[patchI],
                        fields[i].dimensionedInternalField(),
                        fieldDicts.subDict
                        (
                            fieldNames[i]
                        ).subDict("boundaryField").subDict
                        (
                            subMesh().boundary()[patchI].name()
                        )
                    )
                );
            }
        }
    }
}


// Resize map for individual field
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::resizeMap
(
    const label srcIndex,
    const subMeshMapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<subMeshMapper>& boundaryMapper,
    const List<labelListList>& boundaryReverseMaps,
    const List<PtrList<GeomField> >& srcFields,
    GeomField& field
)
{
    // autoMap the internal field
    field.internalField().autoMap(internalMapper);

    // Reverse map for additional cells
    forAll(srcFields, pI)
    {
        // Fetch field for this processor
        const GeomField& srcField = srcFields[pI][srcIndex];

        field.internalField().rmap
        (
            srcField.internalField(),
            internalReverseMaps[pI]
        );
    }

    // Map physical boundary-fields
    forAll(boundaryMapper, patchI)
    {
        // autoMap the patchField
        field.boundaryField()[patchI].autoMap(boundaryMapper[patchI]);

        // Reverse map for additional patch faces
        forAll(srcFields, pI)
        {
            // Fetch field for this processor
            const GeomField& srcField = srcFields[pI][srcIndex];

            field.boundaryField()[patchI].rmap
            (
                srcField.boundaryField()[patchI],
                boundaryReverseMaps[pI][patchI]
            );
        }
    }
}


// Resize all fields in registry
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::resizeMap
(
    const wordList& names,
    const objectRegistry& mesh,
    const subMeshMapper& internalMapper,
    const List<labelList>& internalReverseMaps,
    const PtrList<subMeshMapper>& boundaryMapper,
    const List<labelListList>& boundaryReverseMaps,
    const List<PtrList<GeomField> >& srcFields
)
{
    forAll(names, indexI)
    {
        // Fetch field from registry
        GeomField& field =
        (
            const_cast<GeomField&>
            (
                mesh.lookupObject<GeomField>(names[indexI])
            )
        );

        // Map the field
        coupledInfo<MeshType>::resizeMap
        (
            indexI,
            internalMapper,
            internalReverseMaps,
            boundaryMapper,
            boundaryReverseMaps,
            srcFields,
            field
        );
    }
}


// Resize boundaryFields for all fields in the registry
template <class MeshType>
template <class GeomField>
void coupledInfo<MeshType>::resizeBoundaries
(
    const objectRegistry& mesh,
    const fvBoundaryMesh& boundary
)
{
    typedef typename GeomField::PatchFieldType PatchFieldType;
    typedef typename GeomField::GeometricBoundaryField GeomBoundaryType;

    HashTable<const GeomField*> fields(mesh.lookupClass<GeomField>());

    forAllIter(typename HashTable<const GeomField*>, fields, fIter)
    {
        // Fetch field from registry
        GeomField& field = const_cast<GeomField&>(*fIter());

        GeomBoundaryType& bf = field.boundaryField();

        // Resize boundary
        label nPatches = boundary.size();
        label nOldPatches = field.boundaryField().size();

        // Create a new list of boundaries
        PtrList<PatchFieldType> newbf(nPatches);

        // Existing fields are mapped with new fvBoundaryMesh references
        for (label patchI = 0; patchI < nOldPatches; patchI++)
        {
            label oldPatchSize = bf[patchI].size();

            newbf.set
            (
                patchI,
                PatchFieldType::New
                (
                    bf[patchI],
                    boundary[patchI],
                    field,
                    subMeshMapper(oldPatchSize, identity(oldPatchSize))
                )
            );
        }

        // Size up new patches
        for (label patchI = nOldPatches; patchI < nPatches; patchI++)
        {
            newbf.set
            (
                patchI,
                PatchFieldType::New
                (
                    boundary[patchI].type(),
                    boundary[patchI],
                    field
                )
            );
        }

        // Transfer contents with new patches
        bf.transfer(newbf);
    }
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Disallow default bitwise assignment
template <class MeshType>
void coupledInfo<MeshType>::operator=(const coupledInfo& rhs)
{
    // Check for assignment to self
    if (this == &rhs)
    {
        FatalErrorIn
        (
            "void coupledInfo::operator=(const Foam::coupledInfo&)"
        )
            << "Attempted assignment to self"
            << abort(FatalError);
    }
}


} // End namespace Foam

// ************************************************************************* //
