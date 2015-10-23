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

Class
    Foam::equationSource

Description
    A set of data sources for the equationReader templated by Type.

SourceFiles
    equationSourceI.H
    equationSource.C
    equationSourceTemplates.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationSource_H
#define equationSource_H

#include "word.H"
#include "UPtrList.H"
#include "wordList.H"
#include "PtrList.H"
#include "GeometricFields.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class equationSource Declaration
\*---------------------------------------------------------------------------*/

template <class Type>
class equationSource
{

    // Private data

        //- Name of templated type
        word templateTypeName_;

        //- Individual "Type"s as sources
        UPtrList<const Type> singles_;

        //- Names associated with singles_
        wordList singleNames_;

        //- Dimensions associated with singles_
        PtrList<dimensionSet> singleDimensions_;

        //- Fields of Type - holds the source data for regular
        //  Field<Type>'s, but also can hold source data for GeometricFields
        //      fields_[sourceIndex][geoIndex][cellIndex]
        //  Where:
        //      sourceIndex is the source variable index
        //      geoIndex is to accommodate GeometricFields:
        //          0 = internalField (or standard Field<Type>)
        //          1+ = boundary patch as a Field<Type>
        //      cellIndex is the index number of the Field<Type>
        PtrList<UPtrList<const Field<Type> > > fields_;

        //- Dimensions associated with the fields_
        PtrList<dimensionSet> fieldDimensions_;

        //- Names associated with the fields_
        wordList fieldNames_;

public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct from components
        explicit equationSource
        (
            const word& templateTypeName
        );

    // Destructor
    ~equationSource();

    // Member functions

        // Access

            //- Individual "Type"s as sources
            inline const UPtrList<const Type>& singles() const;

            //- Names associated with singles
            inline const wordList& singleNames() const;

            //- Dimensions associated with singles_
            inline const PtrList<dimensionSet>& singleDimensions() const;

            //- Fields of Type - holds the source data for regular
            inline const PtrList<UPtrList<const Field<Type> > >&
                fields() const;

            //- Dimensions associated with the fields_
            inline const PtrList<dimensionSet>& fieldDimensions() const;

            //- Names associated with the fields_
            inline const wordList& fieldNames() const;

        // Data lookup

            //- True if lookupName is a valid single source
            bool foundSingle(const word& lookupName) const;

            //- True if lookupName is a valid field source
            bool foundField(const word& lookupName) const;

            //- Returns single sourceIndex for lookupName - fails if not found
            label lookupSingle(const word& lookupName) const;

            //- Returns field sourceIndex for lookupName - fails if not found
            label lookupField(const word& lookupName) const;

            //- Returns the number of fields associated with sourceIndex
            label geoSize(label sourceIndex) const;

            //- Returns the field size associated with source and geo indices
            label fieldSize
            (
                label sourceIndex,
                label geoIndex
            ) const;

            //- Return componentIndex for a given component name
            //  Returns -1 if not found (error handled by calling function)
            label lookupComponentIndex(const word& componentName) const;

            //- Returns the number of singles_
            label nSingles() const;

            //- Returns the number of fields_
            label nFields() const;

        // Data retrieval

            //- Retrieve scalar value from singles
            const scalar& singleValue
            (
                label sourceIndex,
                label componentIndex
            ) const;

            //- Retrieve dimensions from singles
            const dimensionSet& singleDimensions(label sourceIndex) const;

            //- Retrieve name associated with singles
            const word& singleName(label sourceIndex) const;

            //- Retrieve scalar value from field
            const scalar& fieldValue
            (
                label sourceIndex,
                label componentIndex,
                label cellIndex,
                label geoIndex
            ) const;

            //- Retrieve component as entire field
            void fullFieldValue
            (
                scalarField& result,
                label sourceIndex,
                label componentIndex,
                label geoIndex
            ) const;

            //- Retrieve dimensions from field
            const dimensionSet& fieldDimensions(label sourceIndex) const;

            //- Retrieve name associated with field
            const word& fieldName(label sourceIndex) const;

        // Adding data sources

            //- Add single source
            void addSource
            (
                const Type& singleIn,
                const word& name,
                dimensionSet dimensions = dimless
            );

            //- Add dimensionedSingle source
            void addSource
            (
                const dimensioned<Type>& dSingleIn
            );

            //- Add field source
            void addSource
            (
                const Field<Type>& fieldIn,
                const word& name,
                dimensionSet dimensions = dimless
            );

            //- Add dimensioned field source
            template<class GeoMesh>
            void addSource
            (
                const DimensionedField<Type, GeoMesh>& dFieldIn
            );

            //- Add geometric field source
            template<template<class> class PatchField, class GeoMesh>
            void addSource
            (
                const GeometricField<Type, PatchField, GeoMesh>& gFieldIn
            );

        // Removing data sources

            //- Remove a single source and reorder index
            void removeSingle(label sourceIndex);

            //- Remove a field source and reorderd index
            void removeField(label sourceIndex);

        // Input / output

            //- Output sources to a dictionary
            dictionary outputDictionary() const;
};

template<>
label equationSource<scalar>::lookupComponentIndex
(
    const word& componentName
) const;

template<>
const scalar& equationSource<scalar>::singleValue
(
    label sourceIndex,
    label componentIndex
) const;

template<>
const scalar& equationSource<scalar>::fieldValue
(
    label sourceIndex,
    label componentIndex,
    label cellIndex,
    label geoIndex
) const;

template<>
void equationSource<scalar>::fullFieldValue
(
    scalarField& result,
    label sourceIndex,
    label componentIndex,
    label geoIndex
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationSourceI.H"
#ifdef NoRepository
#   include "equationSource.C"
//#   include "equationScalarSource.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
