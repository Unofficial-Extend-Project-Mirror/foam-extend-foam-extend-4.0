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

Class
    Foam::equationReader

Description
    The master class for reading equations from a dictionary file.  This class
    holds all the equations, all the data sources, and gives access to all the
    equation parsing / evaluating functions.

SourceFiles
    equationReaderI.H
    equationReaderPointerFunctions.H
    equationReader.C
    equationReaderAssignPointerFunctions.C
    equationReaderCreateMap.C
    equationReaderEvaluate.C
    equationReaderGetSource.C
    equationReaderIO.C
    equationReaderParse.C
    equationReaderTemplates.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationReader_H
#define equationReader_H

#include "dictionary.H"
//#include "dimensionedScalar.H"
#include "UPtrList.H"
#include "equationList.H"
#include "tokenList.H"
#include "labelList.H"
#include "DimensionedField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class equationReader;
// *** Located in equationReaderIO.C ***
Istream& operator>>(Istream&, equationReader&);
Ostream& operator<<(Ostream&, const equationReader&);

/*---------------------------------------------------------------------------*\
                        Class equationReader Declaration
\*---------------------------------------------------------------------------*/

class equationReader
{
    // Private data
    
        //- Dependency backtrace, used for circular-reference detection
        labelList dependents_;
        
        //- Equations
        equationList eqns_;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 2)
        void (Foam::equationReader::*reportOperationFunction_)
        (
            const label&,
            const label&,
            const dimensionedScalar&
        ) const;

        //- Function pointer (used to avoid a conditional at every operation
        //  for debug level 2)
        void (Foam::equationReader::*reportResultFunction_)
        (
            const dimensionedScalar&
        ) const;

        // Sources - these are locations where variables appearing in the
        // equations can be retrieved
        
            //- Variables that can be looked up in a dictionary
            UPtrList<const dictionary> dictSources_;
            
            //- Words to be looked up in a dictionary (equation operation only
            //  gives indices; this is necessary if the source is a dictionary)
            PtrList<word> dictLookups_;
            
            //- Pointers to source scalarLists, used for lists, fields, etc..
            UPtrList<const scalarList> externalScalarLists_;
            
            //- Dimensions associated with the externalScalarLists
            PtrList<dimensionSet> externalScalarListDimensions_;
            
            //- Names associated with the externalScalarLists
            wordList externalScalarListNames_;
            
            //- Current index associated with the externalScalarLists
            labelList externalScalarListIndex_;
            
            //- Pointers to external dimensionedScalars
            UPtrList<const dimensionedScalar> externalDScalars_;
            
            //- Pointers to external scalars
            UPtrList<const scalar> externalScalars_;
            
            //- Names associated with the external scalars
            wordList externalScalarNames_;

            //- Dimensions associated with the external scalars
            PtrList<dimensionSet> externalScalarDimensions_;
            
            //- The output of each equation in the equationList is sent here,
            //  indexed according to the equationList
            UPtrList<scalar> outputScalars_;
            
            //- Names associated with the output scalars
            wordList outputScalarNames_;

            //- Dimensions associated with the output scalars
            PtrList<dimensionSet> outputScalarDimensions_;

            //- Pointers to output scalarLists, used for lists, fields, etc..
            UPtrList<scalarList> outputScalarLists_;
            
            //- Dimensions associated with the outputScalarLists
            PtrList<dimensionSet> outputScalarListDimensions_;
            
            //- Temporary storage during evaluation declared a member
            //  variable as a bug fix for dropped references and pointers
            PtrList<dimensionedScalar> storage_;
            
            //- Internal scalars for storing constants read from equations.
            //  Moved from equation members to equationReader as a bug fix
            PtrList<scalar> internalScalars_;
            
    // Private Member Functions

        //- Disallow default bitwise copy construct
        equationReader(const equationReader&);

        //- Disallow default bitwise assignment
        void operator=(const equationReader&);

    
        // Main parser functions
        
            //- Parse an equation
            //  *** Located in equationReaderParse.C ***
            void parse(label index);

            // Parses a segment of an equation, used after parenthesis 
            // precedence is determined.  Returns the map index that holds the
            // source info for the result
            //  *** Located in equationReaderParse.C ***
            label parseExpression
            (
                label index,
                const tokenList& tl,
                const labelList& opLvl,
                equationOperationList& map,
                const labelList& indices,
                label& storeIndex
            );


        // General parser support functions
            
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
            );

            //- Modify the equation string going in to the parser:
            //  - change ^ into : to allow detection of powers
            //  - add spaces around ( ) + - * / , ^ to prevent expressions from
            //    combining as words
            string stringPreconditioner(const string& rawString);

            // Search through a string and replace all instances of findMe with
            // replaceWith
            void stringReplaceAll
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
                equationOperationList& map,
                labelList& opLvl,
                labelList& pl
            );

            // Returns a labelList with the indices of the maximum operation
            // level.  If more than one grouping exists, returns the first one.
            labelList findMaxOperation
            (
                const labelList& opLvl,
                const labelList& indices
            );

            // Returns a labelList with the indices of the maximum parenthesis
            // level.  If more than one grouping exists, returns the first one.        
            labelList findMaxParenthesis
            (
                const labelList& parenthesisList,
                const labelList& equationIndices
            ) const;

            // Find all negatives, give the upstream token a negative source
            // index, and trim them from the list
            void absorbNegatives
            (
                const label equationIndex,
                const tokenList& tl,
                labelList& eqnIndices,
                labelList& subEqnIndices,
                equationOperationList& map,
                labelList& opLvl
            );
            
            // Sets the first dimensionedScalar equal to the second one without
            // tripping dimension checking errors
            void dsEqual(dimensionedScalar& dso, const dimensionedScalar& dsi);


        // Indexing-related support functions

            //- Given a labelList 'indices', remove all entries from position
            //  'from' to position 'to', but not including 'exceptFor', if
            //  specified.  Search through labelList 'parent', find the same
            //  entries by value (not position), and remove these also.
            void trimListWithParent
            (
                labelList& parent,
                labelList& indices,
                label from,
                label to,
                label exceptFor = -1
            );

            //- Removes indices from a labelList, including from and to
            //  indices, but not including exceptFor index, if specified.
            void trimList
            (
                labelList& indices,
                label from,
                label to,
                label exceptFor = -1
            );

            //- Search through a labelList and return the index with the value
            //  'value'.  Return -1 if failed.
            label findIndex
            (
                const label value,
                const labelList& indexList
            ) const;


        // Data-handling support functions

            //- Return the source list and source index associated with a variable
            //  name.  Searches in this order:
            //  1. externalDScalars;
            //  2. externalScalars; and
            //  3. dictSources.
            //  Returns the source info upon finding it; does not check for
            //  duplicates.
            equationOperation findSource(const word& varName);

            //- Return the dimensionedScalar associated with a source list and
            //  source index
            //  *** Located in equationReaderGetSource.C ***
            dimensionedScalar getSource
            (
                const label equationIndex,
                const label equationOperationIndex,
                const label maxStoreIndex,
                const label storageOffset
            );

            //- Add a constant to the internal scalar source list
            label addInternalScalar(const scalar& value);

            //- Returns true if the stream is a dimensionedScalar dictionary
            //  entry.
            bool isDimensionedScalar(ITstream& is);

            //- Returns true if the stream is a scalar dictionary entry
            bool isScalar(ITstream& is);

            //- Returns true if the stream is an equation entry
            bool isEquation(ITstream& is);


        // Function-pointer functions, used for efficient evaluation
        // equationReader uses dynamically linked function pointers for
        // efficient evaluation and data sourcing.  The prototypes for most of
        // these functions are located in:
        // *** equationReaderPointerFunctions.H ***
#       include "equationReaderPointerFunctions.H"

            //- Link the functions to the pointers
            //  *** Located in equationReaderAssignFunctionPointers.C ***
            void assignFunctionPointers(const label index);

            //- Pointed to by reportOperationFunction_ when debug level is not
            //  high enough - this function does nothing
            void reportOperationDisabled
            (
                const label& index,
                const label& i,
                const dimensionedScalar& ds
            ) const;

            //- Pointed to by reportOperationFunction_ when debug level is
            //  greater than 1.  This function sends operation-by-operation
            //  information to the console
            void reportOperationEnabled
            (
                const label& index,
                const label& i,
                const dimensionedScalar& ds
            ) const;

            //- Pointed to by reportResultFunction_ when debug level is not
            //  high enough - this function does nothing
            void reportResultDisabled
            (
                const dimensionedScalar& ds
            ) const;

            //- Pointed to by reportResultFunction_ when debug level is not
            //  high enough - this function reports the result operation-by-
            //  operation to the console
            void reportResultEnabled
            (
                const dimensionedScalar& ds
            ) const;



public:

    // Static data members
    
//        static const char* const typeName;
    TypeName("equationReader");


    // Constructor (construct null)
    equationReader();

    // Destructor
    virtual ~equationReader();

    // Member functions
    
        // Adding data sources
        
            //- Add a dictionary
            void addDataSource(const dictionary& dict);
            
            //- Add an external scalar
            void addDataSource
            (
                const scalar& value,
                const word& name,
                const dimensionSet& dimensions = dimless
            );
            
            //- Add an external dimensionedScalar
            void addDataSource(const dimensionedScalar& ds);

            //- Add an external scalarList
            void addDataSource
            (
                const scalarList& slist,
                const word& name,
                const dimensionSet& dimensions = dimless
            );

            //- Add an external dimensionedScalarList
            //  *** Located in equationReaderTemplates.C ***
            template<class GeoMesh>
            void addDataSource
            (
                const DimensionedField<scalar, GeoMesh>& dfield
            );
            
            //- Change the index for a scalarList source
            void setListIndex(const word& name, label newIndex);
            
            //- Change the index for all scalarList sources
            void setListIndex(label newIndex);

        // Data source access
        
            //- Equations (const only)
            inline const equationList& eqns() const;
        
            //- Dictionary source pointers
            inline const UPtrList<const dictionary>& dictSources() const;
            
            //- Access to dictionary source pointers
            inline UPtrList<const dictionary>& dictSources();

            //- Dictionary look-up words (const only)
            inline const PtrList<word>& dictLookups() const;

            //- External scalarList pointers
            inline const
                UPtrList<const scalarList>& externalScalarLists() const;

            //- Access to external scalarList pointers
            inline UPtrList<const scalarList>& externalScalarLists();
            
            //- Dimensions associated with external scalarLists
            inline const
                PtrList<dimensionSet>& externalScalarListDimensions() const;

            //- Access to dimensions associated with external scalarLists
            inline PtrList<dimensionSet>& externalScalarListDimensions();
            
            //- Names associated with the external scalarLists
            inline const wordList& externalScalarListNames() const;

            //- Access to names associated with the external scalarLists
            inline wordList& externalScalarListNames();
            
            //- Current index associated with the externalScalarLists
            inline const labelList& externalScalarListIndex() const;

            //- Access to current index associated with the externalScalarLists
            inline labelList& externalScalarListIndex();

            //- External dimensioned scalar data source pointers
            inline const
                UPtrList<const dimensionedScalar>& externalDScalars() const;
            
            //- Access external dimensioned scalar data source pointers
            inline UPtrList<const dimensionedScalar>& externalDScalars();

            //- External scalar data source pointers
            inline const UPtrList<const scalar>& externalScalars() const;
            
            //- Access external scalar data source pointers
            inline UPtrList<const scalar>& externalScalars();

            //- Const access to external scalar names
            inline const wordList& externalScalarNames() const;

            //- Access external scalar names
            inline wordList& externalScalarNames();

            //- Dimensions associated with external scalars
            inline const
                PtrList<dimensionSet>& externalScalarDimensions() const;

            //- Access to dimensions associated with external scalars
            inline PtrList<dimensionSet>& externalScalarDimensions();

            //- Const access to output scalar data pointers
            inline const
                UPtrList<scalar>& outputScalars() const;
            
            //- Access output scalar data pointers
            inline UPtrList<scalar>& outputScalars();

            //- Const access to output scalar names
            inline const
                wordList& outputScalarNames() const;
            
            //- Access output scalar names
            inline wordList& outputScalarNames();

            //- Const access to output scalar dimensions
            inline const
                PtrList<dimensionSet>& outputScalarDimensions() const;
            
            //- Access output scalar dimensions
            inline PtrList<dimensionSet>& outputScalarDimensions();

            //- Const access to pointers to output scalarLists
            inline const UPtrList<scalarList>& outputScalarLists() const;

            //- Access to pointers to output scalarLists
            inline UPtrList<scalarList>& outputScalarLists();
            
            //- Const access to dimensions associated with outputScalarLists
            inline const
                PtrList<dimensionSet>& outputScalarListDimensions() const;

            //- Access to dimensions associated with outputScalarLists
            inline PtrList<dimensionSet>& outputScalarListDimensions();


        // Equations
        
            //- Creates a new equation but does not parse it.  If equation
            //  exists, throws an error.
            void createEquation
            (
                equation eqn
            );

            //- linkOutput functions
            //  These functions assign or reassign an output variable to an
            //  equation - enables the use of update()

                //- Given equation name, and a dimensionedScalar
                void linkOutput
                (
                    const word& eqnName,
                    dimensionedScalar& outputDVar
                );

                //- Given equation name, and dimensionedScalar components
                void linkOutput
                (
                    const word& eqnName,
                    scalar& value,
                    const word& name,
                    const dimensionSet& dimensions = dimless
                );

                //- Given an equation name, and a dimensionedField
                //  *** Located in equationReaderTemplates ***
                template<class GeoMesh>
                void linkOutput
                (
                    const word& eqnName,
                    DimensionedField<scalar, GeoMesh>& dfield
                );

                //- Given an equation name, and scalarList components
                void linkOutput
                (
                    const word& eqnName,
                    scalarList& outputSList,
                    const dimensionSet& dimensions
                );

                //- Given the equation index, and dimensionedScalar
                void linkOutput
                (
                    label index,
                    dimensionedScalar& outputDVar
                );
                
                //- Given the equation index, and dimensionedScalar components
                void linkOutput
                (
                    label index,
                    scalar& value,
                    const word& name,
                    const dimensionSet& dimensions = dimless
                );

                //- Given the equation index, and DimensionedField
                //  *** Located in equationReaderTemplates.C ***
                template<class GeoMesh>
                void linkOutput
                (
                    label index,
                    DimensionedField<scalar, GeoMesh>& dfield
                );

                //- Given the equation index, and scalarList components
                void linkOutput
                (
                    label index,
                    scalarList& outputSList,
                    const dimensionSet& dimensions
                );

            //- readEquation functions - if equation exists, it assumes this is
            //  a reread.  If parameters include output variables, this
            //  function will also perform a linkOutput

                //- Given an equation
                void readEquation
                (
                    equation eqn
                );

                //- Given an equation and linkOutput to a dimensionedScalar
                void readEquation
                (
                    equation eqn,
                    dimensionedScalar& outputDVar
                );

                //- Given an equation and linkOutput to a dimensionedScalar's
                //  components
                void readEquation
                (
                    equation eqn,
                    scalar& value,
                    const word& name,
                    const dimensionSet& dimensions = dimless
                );
                
                //- Given an equation and linkOutput to a DimensionedField
                //  *** Located in equationReaderTemplates.C ***
                template<class GeoMesh>
                void readEquation
                (
                    equation eqn,
                    DimensionedField<scalar, GeoMesh>& dfield
                );

                //- Given an equation and linkOutput to a DimensionedField
                void readEquation
                (
                    equation eqn,
                    scalarList& outputSList,
                    const dimensionSet& dimensions
                );

                //- Read equation from a dictionary
                void readEquation
                (
                    const dictionary& sourceDict,
                    const word& eqnName
                );

                //- Read equation from a dictionary and linkOutput to a
                //  dimensionedScalar
                void readEquation
                (
                    const dictionary& sourceDict,
                    const word& eqnName,
                    dimensionedScalar& outputDVar
                );

                //- Read equation from a dictionary and linkOutput to a
                //  dimensionedScalar's components
                void readEquation
                (
                    const dictionary& sourceDict,
                    const word& eqnName,
                    scalar& value,
                    const word& name,
                    const dimensionSet& dimensions = dimless
                );

            //- Evaluate a single equation, given a variable name, and return
            //  new value.  Do not set the associated output variable
            //  *** Located in equationReaderEvaluate.C ***
            dimensionedScalar evaluate(const word& equationName);
            
            //- Evaluate a single equation, given an index number, and return
            //  new value.  Does not set the associated output variable.
            //  storageOffset is for internal use only, do not use this.
            //  *** Located in equationReaderEvaluate.C ***
            dimensionedScalar evaluate
            (
                const label& index,
                label storageOffset = 0
            );
            
            //- Evaluate the equation and output the result as a Field, cycling
            //  through all the source field indices (setFieldIndex)
            //  Given an index and DimensionedField
            //  *** Located in equationReaderTemplates.C ***
            template<class GeoMesh>
            void evaluateField
            (
                const label& index,
                DimensionedField<scalar, GeoMesh>& dfield
            );
            
            //- Evaluate the equation and output the result as a scalarList,
            //  cycling through all the source field indices (setFieldIndex)
            //  Given an index and scalarList components
            //  *** Located in equationReaderEvaluate.C ***
            void evaluateField
            (
                const label& index,
                scalarList& outputSList,
                dimensionSet& dimensions
            );

            //- Evaluate a single equation, given a variable name, and set the
            //  associated output variable to the new value
            void update(const word& equationName);

            //- Evaluate a single equation, given an index number, and set the
            //  associated output variable to the new value
            void update(const label& index);

            //- Evaluate all equations, setting the associated output variables
            //  as they are calculated
            void update();
            
            //- Returns true if equationName exists in equationList
            bool found(const word& equationName);
            
            //- Returns the index of a given equationName, -1 if not found
            label lookup(const word& equationName);
            
            //- Delete an equation, given a name
            void deleteEquation(const word& equationName);
            
            //- Delete an equation, given an index number
            void deleteEquation(const label& index);
      
        // Input / output related functions
        // *** Located in equationReaderIO.C ***

            //- Output data source information to the Ostream
            Ostream& dataSourceStatus(Ostream& os) const;

            friend Istream& operator>>(Istream&, equationReader&);
            friend Ostream& operator<<(Ostream&, const equationReader&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationReaderI.H"
#include "equationReaderTemplates.C"
#endif
