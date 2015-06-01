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
    Foam::equationReader

Description
    The master class for reading equations from a dictionary file.  This class
    holds all the equations, all the data sources, and gives access to all the
    equation parsing / evaluating functions.

    Has in-depth debug logging to help find problems with equations.  Debug
    switches:

    @verbatim
        0 - None
        1 - Report scalar evaluations (per equation)
        2 - Report scalar evaluations (per operation - very verbose)
        3 - Report dimension evaluations (per equation)
        4 - Report dimension evaluations (per operation - very verbose)
        5 - Combine 1 and 3 above
        6 - Combine 2 and 4 above
    @endverbatim

    Additional index checking implemented in FULLDEBUG mode.

SourceFiles
    Since there are so many files for this object, each function that is not
    implemented in equationReader.C has a comment indicating the associated
    file.

    Inline functions
        equationReaderI.H

    Standard functions
        equationReader.C
        equationReaderAssignFunctionPointers.C
        equationReaderCreateMap.C
        equationReaderEvaluate.C
        equationReaderParse.C

    Templates included with norepository
        equationReaderTemplates.C

    Template specialization
        equationReaders.C

    Input / output functions
        equationReaderIO.C

    Files with functions called through function-pointers
        equationReaderDebugP.H
        equationReaderDebugP.C
        equationReaderEvalDimsP.H
        equationReaderEvalDimsP.C
        equationReaderEvalScalarP.H
        equationReaderEvalScalarP.C
        equationReaderGetSourceDimsP.H
        equationReaderGetSourceDimsP.C
        equationReaderGetSourceScalarP.H
        equationReaderGetSourceScalarP.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationReader_H
#define equationReader_H

#include "versionSpecific.H"
#include "equationReaderVersion.H"
#include "equationSource.H"
#include "equationVariable.H"
#include "dictionary.H"
#include "diagTensor.H"
//#include "dimensionedScalar.H"
#include "UPtrList.H"
//#include "equationList.H"
#include "tokenList.H"
#include "labelList.H"
#include "GeometricFields.H"
//#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class equation;
class equationOperation;

// *** Located in equationReaderIO.C ***
class equationReader;
Istream& operator>>(Istream&, equationReader&);
Ostream& operator<<(Ostream&, const equationReader&);

/*---------------------------------------------------------------------------*\
                        Class equationReader Declaration
\*---------------------------------------------------------------------------*/

class equationReader
//:public PtrList<equation>
{
    // Private data

        //- Equations
        mutable PtrList<equation> eqns_;

        //- Dependency backtrace, used for circular-reference detection
        mutable labelList dependents_;

        // Debug function pointers
        //  equationReader has in-depth debug logging to help users identify
        //  problems with their equations.  The debugging lines are encountered
        //  so frequently that they would seriously impact performance if a
        //  conventional conditional framework is used: "if (debug)".  To
        //  improve execution time, debug function pointers and associated
        //  functions are used.
        // *** Data and prototypes located in equationReaderDebugP.H ***
        // *** Implementation located in equationReaderDebugP.C ***
#       include "equationReaderDebugP.H"

        // Sources - these are locations where variables appearing in the
        // equations can be retrieved

            //- Internal scalars for storing constants read from equations.
            mutable PtrList<scalar> internalScalars_;

            //- Variables that can be looked up in a dictionary
            UPtrList<const dictionary> dictSources_;

            //- Words to be looked up in a dictionary (equation operation only
            //  gives indices; this is necessary if the source is a dictionary)
            mutable PtrList<word> dictLookups_;

            //- Active equation variable sources
            mutable UPtrList<const equationVariable> activeSources_;

            //- List of active equation variable names
            mutable wordList activeSourceNames_;

            //- External scalar sources
            equationSource<scalar> scalarSources_;

            //- External vector sources
            equationSource<vector> vectorSources_;

            //- External tensor sources
            equationSource<tensor> tensorSources_;

            //- External diagTensor sources
            equationSource<diagTensor> diagTensorSources_;

            //- External symmetricTensor sources
            equationSource<symmTensor> symmTensorSources_;

            //- External sphericalTensor sources
            equationSource<sphericalTensor> sphericalTensorSources_;

            //- Current geoIndex of source fields
            mutable label geoIndex_;

            //- Current cellIndex of source fields
            mutable label cellIndex_;

            //- Temporary storage during evaluation
            //mutable fieldSize_;
            mutable scalarField tempSrcField_;
            mutable PtrList<scalarField> storageScalarFields_;
            mutable scalarList storageScalars_;
            mutable PtrList<dimensionSet> storageDims_;

    // Private Member Functions

        //- Disallow default bitwise copy construct
        equationReader(const equationReader&);

        //- Disallow default bitwise assignment
        void operator=(const equationReader&);

        //- Internal dimensions evaluation function
        //  *** Located in equationReaderEvaluate.C ***
        dimensionSet internalEvaluateDimensions
        (
            const label& equationIndex,
            label storageOffset
        ) const;

        //- Internal scalar evaluation function
        //  *** Located in equationReaderEvaluate.C ***
        scalar internalEvaluateScalar
        (
            const label& equationIndex,
            label storageOffset
        ) const;

        //- Internal scalarField evaluation function
        //  *** Located in equationReaderEvaluate.C ***
        void internalEvaluateScalarField
        (
            scalarField& result,
            const label& equationIndex,
            label storageOffset
        ) const;

        //- Internal final dimension checking / error throwing
        //  *** Located in equationReaderEvaluate.C ***
        void checkFinalDimensions
        (
            const label& equationIndex,
            dimensionSet& expectedDimensions,
            const word& outputName
        ) const;


        // Main parser functions

            //- Parse an equation
            //  *** Located in equationReaderParse.C ***
            void parse(label index) const;

            // Parses a segment of an equation, used after parenthesis
            // precedence is determined.  Returns the map index that holds the
            // source info for the result
            //  *** Located in equationReaderParse.C ***
            label parseExpression
            (
                label index,
                const tokenList& tl,
                const labelList& opLvl,
                PtrList<equationOperation>& map,
                const labelList& indices,
                label& storeIndex
            ) const;


        // General parser support functions

            //- Creates a new equation but does not parse it.  If equation
            //  exists, throws an error.  Returns equationIndex.
            label createEquation
            (
                equation& eqn
            ) const;

            //- Post an error related to equation parsing. Reports the exact
            //  position of the error in the equation string.
            void fatalParseError
            (
                const label index,
                const tokenList& tl,
                const label fromToken,
                const label toToken,
                const string& errorIn,
                const OStringStream& description
            ) const;

            //- Modify the equation string going in to the parser:
            //  - change ^ into : to allow detection of powers
            //  - add spaces around ( ) + - * / , ^ to prevent expressions from
            //    combining as words
            static string stringPreconditioner(const string& rawString);

            //- Search through a string and replace all instances of findMe
            //  with replaceWith
            static void stringReplaceAll
            (
                string& working,
                const string& findMe,
                const string& replaceWith
            );

            //- Create a map of: parenthesis levels, operation levels,
            //  functions, and variable sources
            //  *** Located in equationReaderCreateMap.C ***
            void createMap
            (
                const label index,
                const tokenList& tl,
                PtrList<equationOperation>& map,
                labelList& opLvl,
                labelList& pl
            ) const;

            //- Bug fix - once map is created, we find any 'pow(a,b)' instances
            //  and remove the b part to a seperate equation
            void removePowExponents
            (
                const label index,
                tokenList& tl,
                PtrList<equationOperation>& map,
                labelList& opLvl,
                labelList& pl
            ) const;

            //- Returns a labelList with the indices of the maximum operation
            //  level.  If more than one grouping exists, returns the first
            //  one.
            labelList findMaxOperation
            (
                const labelList& opLvl,
                const labelList& indices
            ) const;

            //- Returns a labelList with the indices of the maximum parenthesis
            //  level.  If more than one grouping exists, returns the first
            //  one.
            labelList findMaxParenthesis
            (
                const labelList& parenthesisList,
                const labelList& equationIndices
            ) const;

            //- Find all negatives, give the upstream token a negative source
            //  index, and trim them from the list
            void absorbNegatives
            (
                const label equationIndex,
                const tokenList& tl,
                labelList& eqnIndices,
                labelList& subEqnIndices,
                PtrList<equationOperation>& map,
                labelList& opLvl
            ) const;

            //- Finds the maximum field size that can be returned from an
            //  equation, based on its sources, indexed by geoIndex
            void findMaxFieldSize(const label equationIndex) const;

            //- Sets the first dimensionedScalar equal to the second one
            //  without tripping dimension checking errors
            void dsEqual
            (
                dimensionedScalar& dso,
                const dimensionedScalar& dsi
            ) const;


        // Indexing-related support functions

            //- Given a labelList 'indices', remove all entries from position
            //  'from' to position 'to', but not including 'exceptFor', if
            //  specified.  Search through labelList 'parent', find the same
            //  entries by value (not position), and remove these also.
            static void trimListWithParent
            (
                labelList& parent,
                labelList& indices,
                label from,
                label to,
                label exceptFor = -1
            );

            //- Removes indices from a labelList, including from and to
            //  indices, but not including exceptFor index, if specified.
            static void trimList
            (
                labelList& indices,
                label from,
                label to,
                label exceptFor = -1
            );

            //- Search through a labelList and return the index with the value
            //  'value'.  Return -1 if failed.
            static label findIndex
            (
                const label value,
                const labelList& indexList
            );


        // Data-handling support functions

            //- Return the source list and source index associated with a
            //  variable name.  Searches in this order:
            //  - other known equations
            //  - activeSources;
            //  - scalarSources, then scalarFieldSources;
            //  - vectorSources, then vectorFieldSources;
            //  - tensorSources, then tensorFieldSources;
            //  - diagTensorSources, then diagTensorFieldSources;
            //  - symmTensorSources, then symmTensorFieldSources;
            //  - sphericalTensorSources, then sphericalTensorFieldSources;
            //  - dictSources, which may be either:
            //  -- a scalar;
            //  -- a dimensionedScalar; or
            //  -- a yet unknown equation;
            //  Returns the source info upon finding it; does not check for
            //  duplicates.
            equationOperation findSource(const word& varName) const;

            //- Add a constant to the internal scalar source list
            label addInternalScalar(const scalar& value) const;

        // Function-pointer functions, used for efficient evaluation
        // equationReader uses dynamically linked function pointers for
        // efficient evaluation and data sourcing.  This has the effect of
        // moving all conditionals to the initial parsing stage, allowing the
        // equations to be evaluated without conditonals.  Since conditionals
        // are slow, execution speed approaches that of hard-coded equations.
        //
        // Data sourcing functions for dimensionSet evaluation:
        // *** Prototypes located in equationReaderGetSourceDimsP.H ***
        // *** Implemented in equationReaderGetSourceDimsP.C ***
        //
        // Data sourcing functions for scalar evaluation:
        // *** Prototypes located in equationReaderGetSourceScalarP.H ***
        // *** Implemented in equationReaderGetSourceScalarP.C ***
        //
        // Data sourcing functions for scalarField evaluation:
        // *** Prototypes located in equationReaderGetSourceScalarFieldP.H ***
        // *** Implemented in equationReaderGetSourceScalarFieldP.C ***
        //
        // Evaluation functions for dimensionSet evaluation:
        // *** Prototypes located in equationReaderEvalDimsP.H ***
        // *** Implemented in equationReaderEvalDimsP.C ***
        //
        // Evaluation functions for scalar evaluation:
        // *** Prototypes located in equationReaderEvalScalarP.H ***
        // *** Implemented in equationReaderEvalScalarP.C ***
        //
        // Evaluation functions for scalarField evaluation:
        // *** Prototypes located in equationReaderEvalScalarFieldP.H ***
        // *** Implemented in equationReaderEvalScalarFieldP.C ***
#       include "equationReaderGetSourceDimsP.H"
#       include "equationReaderGetSourceScalarP.H"
#       include "equationReaderGetSourceScalarFieldP.H"
#       include "equationReaderEvalDimsP.H"
#       include "equationReaderEvalScalarP.H"
#       include "equationReaderEvalScalarFieldP.H"

            //- Link the functions to the pointers
            //  *** Located in equationReaderAssignFunctionPointers.C ***
            void assignFunctionPointers(const label index) const;

        // Functions to emulate "is a" PtrList<equation>, rather than "has a"
        // (because the equations have to be mutable)

            //- Set size of equation list (private)
            void setSize(const label newSize) const;

            //- Set the an element of the equation list (private)
            void set(const label index, equation * eqn) const;


public:

    // Static data members

    TypeName("equationReader");


    // Constructor (construct null)
    equationReader(bool showSplash = true);

    // Destructor
    virtual ~equationReader();

    // Member functions

        // Access

            //- Version number
            word version() const;

            //- Equation list
            inline const PtrList<equation>& eqns() const;

            //- Internal scalars
            //  *** Located in equationReaderI.H ***
            inline const PtrList<scalar>& internalScalars() const;

            //- Dictionary sources
            //  *** Located in equationReaderI.H ***
            inline const UPtrList<const dictionary>& dictSources() const;

            //- Dictionary lookup words
            //  *** Located in equationReaderI.H ***
            inline const PtrList<word>& dictLookups() const;

            //- Active equation variable sources
            //  *** Located in equationReaderI.H ***
            inline const UPtrList<const equationVariable>&
                activeSources() const;

            //- Active equation variable source names
            //  *** Located in equationReaderI.H ***
            inline const wordList& activeSourceNames() const;

            //- Access scalar sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<scalar>& scalarSources() const;
            inline equationSource<scalar>& scalarSources();

            //- Access vector sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<vector>& vectorSources() const;
            inline equationSource<vector>& vectorSources();

            //- Access tensor sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<tensor>& tensorSources() const;
            inline equationSource<tensor>& tensorSources();

            //- Access diagTensor sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<diagTensor>& diagTensorSources() const;
            inline equationSource<diagTensor>& diagTensorSources();

            //- Access symmetricTensor sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<symmTensor>& symmTensorSources() const;
            inline equationSource<symmTensor>& symmTensorSources();

            //- Access sphericalTensor sources
            //  *** Located in equationReaderI.H ***
            inline const equationSource<sphericalTensor>&
                sphericalTensorSources() const;
            inline equationSource<sphericalTensor>& sphericalTensorSources();

            //- Current geoIndex with the fields
            //  *** Located in equationReaderI.H ***
            inline const label& geoIndex() const;

            //- Change the index for fields
            //  *** Located in equationReaderI.H ***
            inline void setGeoIndex(label newIndex);

            //- Current cellIndex with the fields
            //  *** Located in equationReaderI.H ***
            inline const label& cellIndex() const;

            //- Change the index for fields
            //  *** Located in equationReaderI.H ***
            inline void setCellIndex(label newIndex);

        // Data lookup

            //- Returns true if equationName exists in equation list
            bool found(const word& equationName);

            //- Returns the index of a given equationName, -1 if not found
            label lookup(const word& equationName) const;

        // Adding data sources

            //- Add a dictionary
            void addSource(const dictionary& dict);

            //- Add an active variable
            void addSource(const equationVariable& aVar);

        // Equations

            //- Read an equation, given the equation.
            //  okayToReread: false=throw error if equation already exists
            //                true=reread the equation
            //  Returns equationIndex that is assigned
            label readEquation
            (
                equation eqn,
                bool okayToReread = false
            );

            //- Read an equation from a dictionary
            //  okayToReread: false=throw error if equation already exists
            //                true=reread the equation
            //  Returns equationIndex that is assigned
            label readEquation
            (
                const dictionary& sourceDict,
                const word& eqnName,
                bool okayToReread = false
            );

            //- Force an equation to re-parse
            void clearEquation(const label equationIndex) const;

            //- Delete an equation, given a name
            void deleteEquation(const word& equationName);

            //- Delete an equation, given an index number
            void deleteEquation(const label& index);

        // Evaluate equations

            //- Return a scalar, given the equation name
            //  *** Located in equationReaderEvaluate.C ***
            scalar evaluateScalar
            (
                const word& equationName,
                const label cellIndex = 0,
                const label geoIndex = 0
            ) const;

            //- Return a scalar, given the equation index
            //  *** Located in equationReaderEvaluate.C ***
            scalar evaluateScalar
            (
                const label equationIndex,
                const label cellIndex = 0,
                const label geoIndex = 0
            ) const;

            //- Return the resulting dimensions, given an equation name
            //  *** Located in equationReaderEvaluate.C ***
            dimensionSet evaluateDimensions(const word& equationName) const;

            //- Return the resulting dimensions, given an equation index
            //  *** Located in equationReaderEvaluate.C ***
            dimensionSet evaluateDimensions(const label equationIndex) const;

            //- Return a dimensionedScalar, given the equation name
            //  *** Located in equationReaderEvaluate.C ***
            dimensionedScalar evaluateDimensionedScalar
            (
                const word& equationName,
                const label cellIndex = 0,
                const label geoIndex = 0
            ) const;

            //- Return a dimensionedScalar, given the equation index
            //  *** Located in equationReaderEvaluate.C ***
            dimensionedScalar evaluateDimensionedScalar
            (
                const label equationIndex,
                const label cellIndex = 0,
                const label geoIndex = 0
            ) const;

            //- Fill a given scalarField, given an equation name
            //  *** Located in equationReaderEvaluate.C ***
            void evaluateScalarField
            (
                scalarField& resultField,
                const word& equationName,
                const label geoIndex = 0
            ) const;

            //- Fill a given scalarField, given an equation index
            //  *** Located in equationReaderEvaluate.C ***
            void evaluateScalarField
            (
                scalarField& resultField,
                const label equationIndex,
                const label geoIndex = 0
            ) const;

            //- Fill the named component of a given Field<Type>, given an
            //  equation name
            //  *** Located in equationReaderTemplates.C ***
            template<class Type>
            void evaluateTypeField
            (
                Field<Type>& resultField,
                const word& componentName,
                const word& equationName,
                const label geoIndex = 0
            ) const;

            //- Fill the given component of a given Field<Type>, given an
            //  equation index
            //  *** Located in equationReaderTemplates.C ***
            template<class Type>
            void evaluateTypeField
            (
                Field<Type>& resultField,
                const label componentIndex,
                const label equationIndex,
                const label geoIndex = 0
            ) const;

            //- Fill a given dimensioned scalarField, given an equation name
            //  *** Located in equationReaderTemplates.C ***
            template<class GeoMesh>
            void evaluateDimensionedScalarField
            (
                DimensionedField<scalar, GeoMesh>& resultDField,
                const word& equationName,
                const label geoIndex = 0
            ) const;

            //- Fill a given dimensioned scalarField, given an equation index
            //  *** Located in equationReaderTemplates.C ***
            template<class GeoMesh>
            void evaluateDimensionedScalarField
            (
                DimensionedField<scalar, GeoMesh>& resultDField,
                const label equationIndex,
                const label geoIndex = 0
            ) const;

            //- Fill the named component of a given dimensionedField<Type>,
            //  given an equation name
            //  *** Located in equationReaderTemplates.C ***
            template<class Type, class GeoMesh>
            void evaluateDimensionedTypeField
            (
                DimensionedField<Type, GeoMesh>& resultDField,
                const word& componentName,
                const word& equationName,
                const label geoIndex = 0
            ) const;

            //- Fill the given component of a given dimensionedField<Type>,
            //  given an equation index
            //  *** Located in equationReaderTemplates.C ***
            template<class Type, class GeoMesh>
            void evaluateDimensionedTypeField
            (
                DimensionedField<Type, GeoMesh>& resultDField,
                const label componentIndex,
                const label equationIndex,
                const label geoIndex = 0
            ) const;

            //- Fill a given GeometricScalarField, given an equation name
            //  *** Located in equationReaderTemplates.C ***
            template<template<class> class PatchField, class GeoMesh>
            void evaluateGeometricScalarField
            (
                GeometricField<scalar, PatchField, GeoMesh>& resultGField,
                const word& equationName
            ) const;

            //- Fill a given GeometricField, given an equation index
            //  *** Located in equationReaderTemplates.C ***
            template <template<class> class PatchField, class GeoMesh>
            void evaluateGeometricScalarField
            (
                GeometricField<scalar, PatchField, GeoMesh>& resultGField,
                const label equationIndex
            ) const;

            //- Fill the named component of a given GeometricField<Type>,
            //  given an equation name
            //  *** Located in equationReaderTemplates.C ***
            template
            <
                class Type, template<class> class PatchField, class GeoMesh
            >
            void evaluateGeometricTypeField
            (
                GeometricField<Type, PatchField, GeoMesh>& resultGField,
                const word& componentName,
                const word& equationName
            ) const;

            //- Fill the given component of a given GeometricField<Type>,
            //  given an equation index
            //  *** Located in equationReaderTemplates.C ***
            template
            <
                class Type, template<class> class PatchField, class GeoMesh
            >
            void evaluateGeometricTypeField
            (
                GeometricField<Type, PatchField, GeoMesh>& resultGField,
                const label componentIndex,
                const label equationIndex
            ) const;

        // Tools

            //- Returns true if the stream is a dimensionedScalar dictionary
            //  entry.
            static bool isDimensionedScalar(ITstream& is);

            //- Returns true if the stream is a word dictionary entry
            static bool isWord(ITstream& is);

            //- Returns true if the stream is a scalar dictionary entry
            static bool isScalar(ITstream& is);

            //- Returns true if the stream is a nameless dimensionedScalar
            //  dictionary entry
            static bool isNamelessDimensionedScalar(ITstream& is);

            //- Returns true if the stream is an equation entry
            static bool isEquation(ITstream& is);

        // Functions to emulate "is a" PtrList<equation>, rather than "has a"
        // (because the equations have to be mutable)

            //- Index operator
            const equation& operator[](const label equationIndex) const;

            //- Index operator
            equation& operator[](const label equationIndex);

            //- Return size of equation list
            label size() const;

        // Input / output related functions
        // *** Located in equationReaderIO.C ***

            //- Output data source information to the Ostream
            Ostream& dataSourceStatus(Ostream& os) const;

            friend Istream& operator>>(Istream&, equationReader&);
            friend Ostream& operator<<(Ostream&, const equationReader&);
};

// Scalar specialisation
template<>
void equationReader::evaluateTypeField
(
    scalarField& resultField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const;

// Scalar specialisation
template<>
void equationReader::evaluateTypeField
(
    scalarField& resultField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const;

/*
// Scalar specialisation
template<class GeoMesh>
void evaluateDimensionedTypeField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const word& componentName,
    const word& equationName,
    const label geoIndex
) const;

// Scalar specialisation
template<class GeoMesh>
void evaluateDimensionedTypeField
(
    DimensionedField<scalar, GeoMesh>& resultDField,
    const label componentIndex,
    const label equationIndex,
    const label geoIndex
) const;

// Scalar specialisation
template<template<class> class PatchField, class GeoMesh>
evaluateGeometricTypeField
(
    GeometricField
    <
        scalar, template<class> class PatchField, class GeoMesh
    >& resultGField,
    const word& componentName,
    const word& equationName
) const;

// Scalar specialisation
template<template<class> class PatchField, class GeoMesh>
evaluateGeometricTypeField
(
    GeometricField
    <
        scalar, template<class> class PatchField, class GeoMesh
    >& resultGField,
    const label componentIndex,
    const label equationIndex
) const;
*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationReaderI.H"
#ifdef NoRepository
#   include "equationReaderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
