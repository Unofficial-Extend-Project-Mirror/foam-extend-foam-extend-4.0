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

#include "Time.H"
#include "coupledInfo.H"
#include "dynamicTopoFvMesh.H"
#include "emptyFvPatchFields.H"
#include "emptyFvsPatchFields.H"
#include "fixedValueFvPatchFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Constructor for coupledInfo
coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
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


coupledInfo::coupledInfo
(
    const dynamicTopoFvMesh& mesh,
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
coupledInfo::subMeshMapper::subMeshMapper
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

const dynamicTopoFvMesh& coupledInfo::baseMesh() const
{
    return mesh_;
}


void coupledInfo::setMesh
(
    label index,
    dynamicTopoFvMesh* mesh
)
{
    subMesh_.set(mesh);
}


dynamicTopoFvMesh& coupledInfo::subMesh()
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("dynamicTopoFvMesh& coupledInfo::subMesh()")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


const dynamicTopoFvMesh& coupledInfo::subMesh() const
{
    if (!subMesh_.valid())
    {
        FatalErrorIn("const dynamicTopoFvMesh& coupledInfo::subMesh() const")
            << " Sub-mesh pointer has not been set."
            << abort(FatalError);
    }

    return subMesh_();
}


bool coupledInfo::builtMaps() const
{
    return builtMaps_;
}


void coupledInfo::setBuiltMaps()
{
    builtMaps_ = true;
}


coupleMap& coupledInfo::map()
{
    return map_;
}


const coupleMap& coupledInfo::map() const
{
    return map_;
}


label coupledInfo::masterFaceZone() const
{
    return masterFaceZone_;
}


label coupledInfo::slaveFaceZone() const
{
    return slaveFaceZone_;
}


// Set subMesh centres
void coupledInfo::setCentres(PtrList<volVectorField>& centres) const
{
    // Fetch reference to subMesh
    const dynamicTopoFvMesh& mesh = subMesh();

    // Set size
    centres.setSize(1);

    vectorField Cv(mesh.cellCentres());
    vectorField Cf(mesh.faceCentres());

    // Create and map the patch field values
    label nPatches = mesh.boundary().size();

    // Create field parts
    PtrList<fvPatchField<vector> > volCentrePatches(nPatches);

    // Over-ride and set all patches to fixedValue
    for (label patchI = 0; patchI < nPatches; patchI++)
    {
        volCentrePatches.set
        (
            patchI,
            new fixedValueFvPatchField<vector>
            (
                mesh.boundary()[patchI],
                DimensionedField<vector, volMesh>::null()
            )
        );

        // Slice field to patch (forced assignment)
        volCentrePatches[patchI] ==
        (
            mesh.boundaryMesh()[patchI].patchSlice(Cf)
        );
    }

    // Set the cell-centres pointer.
    centres.set
    (
        0,
        new volVectorField
        (
            IOobject
            (
                "cellCentres",
                mesh.time().timeName(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh,
            dimLength,
            SubField<vector>(Cv, mesh.nCells()),
            volCentrePatches
        )
    );
}


// Subset volume field
template <class Type>
tmp<GeometricField<Type, fvPatchField, volMesh> >
coupledInfo::subSetVolField
(
    const GeometricField<Type, fvPatchField, volMesh>& fld
) const
{
    // Create and map the internal-field values
    Field<Type> internalField
    (
        fld.internalField(),
        map().cellMap()
    );

    // Create and map the patch field values
    label nPatches = subMesh().boundary().size();
    PtrList<fvPatchField<Type> > patchFields(nPatches);

    forAll(patchFields, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            patchFields.set
            (
                patchI,
                new emptyFvPatchField<Type>
                (
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, volMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                fvPatchField<Type>::New
                (
                    fld.boundaryField()[patchI],
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, volMesh>::null(),
                    subMeshMapper(*this, patchI)
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeometricField<Type, fvPatchField, volMesh> > subFld
    (
        new GeometricField<Type, fvPatchField, volMesh>
        (
            IOobject
            (
                "subField_" + fld.name(),
                subMesh().time().timeName(),
                subMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh(),
            fld.dimensions(),
            internalField,
            patchFields
        )
    );

    return subFld;
}


// Subset surface field
template <class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
coupledInfo::subSetSurfaceField
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& fld
) const
{
    // Create and map the internal-field values
    Field<Type> internalField
    (
        fld.internalField(),
        SubList<label>
        (
            map().faceMap(),
            subMesh().nInternalFaces()
        )
    );

    // Create and map the patch field values
    label nPatches = subMesh().boundary().size();
    PtrList<fvsPatchField<Type> > patchFields(nPatches);

    forAll(patchFields, patchI)
    {
        if (patchI == (nPatches - 1))
        {
            // Artificially set last patch
            patchFields.set
            (
                patchI,
                new emptyFvsPatchField<Type>
                (
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null()
                )
            );
        }
        else
        {
            patchFields.set
            (
                patchI,
                fvsPatchField<Type>::New
                (
                    fld.boundaryField()[patchI],
                    subMesh().boundary()[patchI],
                    DimensionedField<Type, surfaceMesh>::null(),
                    subMeshMapper(*this, patchI)
                )
            );
        }
    }

    // Create new field from pieces
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > subFld
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "subField_" + fld.name(),
                subMesh().time().timeName(),
                subMesh(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            subMesh(),
            fld.dimensions(),
            internalField,
            patchFields
        )
    );

    return subFld;
}


template <class Type>
void coupledInfo::mapVolField
(
    const wordList& fieldNames,
    const word& fieldType,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        const GeometricField<Type, fvPatchField, volMesh>& fld =
        (
            mesh_.lookupObject
            <
                GeometricField<Type, fvPatchField, volMesh>
            >(fieldNames[i])
        );

        tmp<GeometricField<Type, fvPatchField, volMesh> > tsubFld =
        (
            subSetVolField(fld)
        );

        // Send field through stream
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


template <class Type>
void coupledInfo::mapSurfaceField
(
    const wordList& fieldNames,
    const word& fieldType,
    OSstream& strStream
) const
{
    strStream
        << fieldType << token::NL
        << token::BEGIN_BLOCK << token::NL;

    forAll(fieldNames, i)
    {
        const GeometricField<Type, fvsPatchField, surfaceMesh>& fld =
        (
            mesh_.lookupObject
            <
                GeometricField<Type, fvsPatchField, surfaceMesh>
            >(fieldNames[i])
        );

        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsubFld =
        (
            subSetSurfaceField(fld)
        );

        // Send field through stream
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


// Set volume field pointer from input dictionary
template <class GeomField>
void coupledInfo::setField
(
    const wordList& fieldNames,
    const dictionary& fieldDicts,
    PtrList<GeomField>& fields
) const
{
    // Size up the pointer list
    fields.setSize(fieldNames.size());

    forAll(fieldNames, i)
    {
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
                fieldDicts.subDict(fieldNames[i])
            )
        );
    }
}


template <class GeomField>
void coupledInfo::resizeMap
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
template <class GeomField>
void coupledInfo::resizeMap
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
        coupledInfo::resizeMap
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


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void coupledInfo::operator=(const coupledInfo& rhs)
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
