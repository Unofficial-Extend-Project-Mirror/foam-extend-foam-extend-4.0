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

#include "equationSource.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::equationSource<Type>::equationSource
(
    const word& templateTypeName
)
:
    templateTypeName_(templateTypeName)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::equationSource<Type>::~equationSource()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
bool Foam::equationSource<Type>::foundSingle(const word& lookupName) const
{
    forAll(singleNames_, sourceIndex)
    {
        if (singleNames_[sourceIndex] == lookupName)
        {
            return true;
        }
    }
    return false;
}


template<class Type>
bool Foam::equationSource<Type>::foundField(const word& lookupName) const
{
    forAll(fieldNames_, sourceIndex)
    {
        if (fieldNames_[sourceIndex] == lookupName)
        {
            return true;
        }
    }
    return false;
}


template<class Type>
Foam::label Foam::equationSource<Type>::lookupSingle
(
    const word& lookupName
) const
{
    forAll(singleNames_, sourceIndex)
    {
        if (singleNames_[sourceIndex] == lookupName)
        {
            return sourceIndex;
        }
    }
    FatalErrorIn("equationSource::lookupSingle")
        << lookupName << " is not a valid " << templateTypeName_ << " single "
        << "data source.  Valid names are:" << token::NL << singleNames_
        << abort(FatalError);
    return -1;
}


template<class Type>
Foam::label Foam::equationSource<Type>::lookupField
(
    const word& lookupName
) const
{
    forAll(fieldNames_, sourceIndex)
    {
        if (fieldNames_[sourceIndex] == lookupName)
        {
            return sourceIndex;
        }
    }
    FatalErrorIn("equationSource::lookupField")
        << lookupName << " is not a valid " << templateTypeName_ << "Field "
        << "data source.  Valid names are:" << token::NL << fieldNames_
        << abort(FatalError);
    return -1;
}


template<class Type>
Foam::label Foam::equationSource<Type>::geoSize(label sourceIndex) const
{
    return fields_[sourceIndex].size();
}


template<class Type>
Foam::label Foam::equationSource<Type>::fieldSize
(
    label sourceIndex,
    label geoIndex
) const
{
    return fields_[sourceIndex][geoIndex].size();
}


template<class Type>
Foam::label Foam::equationSource<Type>::lookupComponentIndex
(
    const word& componentName
) const
{
    // Because Type.size() is not static
    Type dummy;

    for (label compIndex(0); compIndex < dummy.size(); compIndex++)
    {
        if (Type::componentNames[compIndex] == componentName)
        {
            return compIndex;
        }
    }
    return -1;
}


template<class Type>
Foam::label Foam::equationSource<Type>::nSingles() const
{
    return singles_.size();
}


template<class Type>
Foam::label Foam::equationSource<Type>::nFields() const
{
    return fields_.size();
}


template<class Type>
const Foam::scalar& Foam::equationSource<Type>::singleValue
(
    label sourceIndex,
    label componentIndex
) const
{
    return singles_[sourceIndex][componentIndex];
}


template<class Type>
const Foam::dimensionSet& Foam::equationSource<Type>::singleDimensions
(
    label sourceIndex
) const
{
    return singleDimensions_[sourceIndex];
}


template<class Type>
const Foam::word& Foam::equationSource<Type>::singleName
(
    label sourceIndex
) const
{
    return singleNames_[sourceIndex];
}


template<class Type>
const Foam::scalar& Foam::equationSource<Type>::fieldValue
(
    label sourceIndex,
    label componentIndex,
    label cellIndex,
    label geoIndex
) const
{
    return fields_[sourceIndex][geoIndex][cellIndex][componentIndex];
}


template<class Type>
void Foam::equationSource<Type>::fullFieldValue
(
    scalarField& result,
    label sourceIndex,
    label componentIndex,
    label geoIndex
) const
{
    const Field<Type>& fieldRef(fields_[sourceIndex][geoIndex]);
    //result.setSize(fieldRef.size());
    forAll(result, cellIndex)
    {
        result[cellIndex] = fieldRef[cellIndex][componentIndex];
    }
}


template<class Type>
const Foam::dimensionSet& Foam::equationSource<Type>::fieldDimensions
(
    label sourceIndex
) const
{
    return fieldDimensions_[sourceIndex];
}


template<class Type>
const Foam::word& Foam::equationSource<Type>::fieldName
(
    label sourceIndex
) const
{
    return fieldNames_[sourceIndex];
}


template<class Type>
void Foam::equationSource<Type>::addSource
(
    const Type& singleIn,
    const word& name,
    dimensionSet dimensions
)
{
    label newIndex(singles_.size());
    singles_.setSize(newIndex + 1);
    singleDimensions_.setSize(newIndex + 1);
    singleNames_.setSize(newIndex + 1);

    singles_.set
    (
        newIndex,
        &singleIn
    );
    singleDimensions_.set
    (
        newIndex,
        new dimensionSet(dimensions)
    );
    singleNames_[newIndex] = name;
}


template<class Type>
void Foam::equationSource<Type>::addSource
(
    const dimensioned<Type>& dSingleIn
)
{
    addSource(dSingleIn.value(), dSingleIn.name(), dSingleIn.dimensions());
}


template<class Type>
void Foam::equationSource<Type>::addSource
(
    const Field<Type>& fieldIn,
    const word& name,
    dimensionSet dimensions
)
{
    label newIndex(fields_.size());
    fields_.setSize(newIndex + 1);
    fieldDimensions_.setSize(newIndex + 1);
    fieldNames_.setSize(newIndex + 1);

    fields_.set
    (
        newIndex,
        new UPtrList<const Field<Type> >(1)
    );
    fields_[newIndex].set
    (
        0,
        &fieldIn
    );
    fieldDimensions_.set
    (
        newIndex,
        new dimensionSet(dimensions)
    );
    fieldNames_[newIndex] = name;
}


template<class Type>
template<class GeoMesh>
void Foam::equationSource<Type>::addSource
(
    const DimensionedField<Type, GeoMesh>& dFieldIn
)
{
    addSource
    (
        dFieldIn.field(),
        dFieldIn.name(),
        dFieldIn.dimensions()
    );
}


template<class Type>
template<template<class> class PatchField, class GeoMesh>
void Foam::equationSource<Type>::addSource
(
    const GeometricField<Type, PatchField, GeoMesh>& gFieldIn
)
{
    label newIndex(fields_.size());
    label newGeoIndex(gFieldIn.boundaryField().size() + 1);

    fields_.setSize(newIndex + 1);
    fieldDimensions_.setSize(newIndex + 1);
    fieldNames_.setSize(newIndex + 1);

    // Set dimensions
    fieldDimensions_.set
    (
        newIndex,
        new dimensionSet(gFieldIn.dimensions())
    );

    // Set name
    fieldNames_[newIndex] = gFieldIn.name();

    // Create fields pointer object
    fields_.set
    (
        newIndex,
        new UPtrList<const Field<Type> >(newGeoIndex)
    );

    // Set internal field
    fields_[newIndex].set
    (
        0,
        &gFieldIn.internalField()
    );

    // Set boundary fields
    forAll(gFieldIn.boundaryField(), patchI)
    {
        fields_[newIndex].set
        (
            patchI + 1,
            &gFieldIn.boundaryField()[patchI]
        );
    }
}


template<class Type>
void Foam::equationSource<Type>::removeSingle(label sourceIndex)
{
#   ifdef FULLDEBUG
        if ((sourceIndex < 0) || (sourceIndex >= nSingles()))
        {
            FatalErrorIn("equationSource::removeSingle")
                << "sourceIndex out of range (0.." << nSingles() << ")"
                << abort(FatalError);
        }
#   endif

    for (label i(sourceIndex); i < (singles_.size() - 1); i++)
    {
        singles_[i] = singles_[i + 1];
        singleDimensions_[i] = singleDimensions_[i + 1];
        singleNames_[i] = singleNames_[i + 1];
    }

    /*
    labelList oldToNew(singles_.size());
    for (label i(0); i < sourceIndex; i++)
    {
        oldToNew[i] = i;
    }
    for (label i(sourceIndex); i < (singles_.size() - 1); i++)
    {
        oldToNew[i] = i + 1;
    }
    oldToNew[singles_.size() - 1] = sourceIndex;

    singles_.reorder(oldToNew);
    singleDimensions_.reorder(oldToNew);
    singleNames_.reorder(oldToNew);*/

    label newSize(singles_.size() - 1);
    singles_.setSize(newSize);
    singleDimensions_.setSize(newSize);
    singleNames_.setSize(newSize);
}


template<class Type>
void Foam::equationSource<Type>::removeField(label sourceIndex)
{
#   ifdef FULLDEBUG
        if ((sourceIndex < 0) || (sourceIndex >= nFields()))
        {
            FatalErrorIn("equationSource::removeSingle")
                << "sourceIndex out of range (0.." << nFields() << ")"
                << abort(FatalError);
        }
#   endif

    for (label i(sourceIndex); i < (fields_.size() - 1); i++)
    {
        fields_[i] = fields_[i + 1];
        fieldDimensions_[i] = fieldDimensions_[i + 1];
        fieldNames_[i] = fieldNames_[i + 1];
    }

    /*
    labelList oldToNew(fields_.size());
    for (label i(0); i < sourceIndex; i++)
    {
        oldToNew[i] = i;
    }
    for (label i(sourceIndex); i < (fields_.size() - 1); i++)
    {
        oldToNew[i] = i + 1;
    }
    oldToNew[fields_.size() - 1] = sourceIndex;

    fields_.reorder(oldToNew);
    fieldDimensions_.reorder(oldToNew);
    fieldNames_.reorder(oldToNew);
    */
    label newSize(fields_.size() - 1);
    fields_.setSize(newSize);
    fieldDimensions_.setSize(newSize);
    fieldNames_.setSize(newSize);
}


template<class Type>
Foam::dictionary Foam::equationSource<Type>::outputDictionary() const
{
    dictionary returnMe;
    dictionary singlesDict;
    forAll(singles_, sourceIndex)
    {
        singlesDict.set
        (
            singleNames_[sourceIndex],
            dimensioned<Type>
            (
                singleNames_[sourceIndex],
                singleDimensions_[sourceIndex],
                singles_[sourceIndex]
            )
        );
    }
    returnMe.set
    (
        keyType(word("single(" + templateTypeName_ + "s)")),
        singlesDict
    );

    dictionary fieldsDict;
    forAll(fields_, sourceIndex)
    {
        dictionary tempDict;
        tempDict.set("dimensions", fieldDimensions_[sourceIndex]);
        forAll(fields_[sourceIndex], geoIndex)
        {
            tempDict.set
            (
                word("field" + name(geoIndex)),
                fields_[sourceIndex][geoIndex]
            );
        }
        fieldsDict.set(fieldNames_[sourceIndex], tempDict);
    }
    returnMe.set
    (
        keyType(word("field(" + templateTypeName_ + "s)")),
        fieldsDict
    );

    return returnMe;
}

// ************************************************************************* //
