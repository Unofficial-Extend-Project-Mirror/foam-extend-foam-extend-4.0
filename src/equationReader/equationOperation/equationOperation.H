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
    Foam::equationOperation

Description
    Defines a single operation to be performed in sequence while evaluating an
    equation read from a dictionary.

SourceFiles
    equationOperationI.H
    equationOperation.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equationOperation_H
#define equationOperation_H

#include "word.H"
#include "label.H"
#include "IOstreams.H"
#include "dimensionedScalar.H"
// #include "Istream.H"
// #include "Ostream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class equationOperation;
class equationReader;

// Friend Operators

int operator==(const equationOperation&, const equationOperation&);
int operator!=(const equationOperation&, const equationOperation&);

// Friend IOstream Operators

Istream& operator>>(Istream&, equationOperation&);
Ostream& operator<<(Ostream&, const equationOperation&);



/*---------------------------------------------------------------------------*\
                        Class equationOperation Declaration
\*---------------------------------------------------------------------------*/

class equationOperation
{

public:

        enum sourceTypeEnum
        {
            stnone,
            ststorage,
            stactiveSource,
            stequation,
            stinternalScalar,
            stdictSource,
            stscalarSource,
            stscalarFieldSource,
            stvectorSource,
            stvectorFieldSource,
            sttensorSource,
            sttensorFieldSource,
            stdiagTensorSource,
            stdiagTensorFieldSource,
            stsymmTensorSource,
            stsymmTensorFieldSource,
            stsphericalTensorSource,
            stsphericalTensorFieldSource
        };

        enum operationType
        {
            otnone,
            otretrieve,
            otstore,
            otplus,
            otminus,
            ottimes,
            otdivide,
            otpow,
            otsign,
            otpos,
            otneg,
            otmag,
            otlimit,
            otminMod,
            otsqrtSumSqr,
            otsqr,
            otpow3,
            otpow4,
            otpow5,
            otpow6,
            otinv,
            otsqrt,
            otcbrt,
            othypot,
            otexp,
            otlog,
            otlog10,
            otsin,
            otcos,
            ottan,
            otasin,
            otacos,
            otatan,
            otatan2,
            otsinh,
            otcosh,
            ottanh,
            otasinh,
            otacosh,
            otatanh,
            oterf,
            oterfc,
            otlgamma,
            otj0,
            otj1,
            otjn,
            oty0,
            oty1,
            otyn,
            otmax,
            otmin,
            otstabilise
        };

private:

    // Source to read the data from
    sourceTypeEnum source_;

    // Index in the field where the data is located.  Note, the
    // equationOperation lists created by equationReader make this variable a
    // 1-indexed (i.e. starts from 1, not zero) in order to use its sign to
    // store the sign of the variable.  The sourceTypes are zero-indexed, so
    // equationReader will constantly be adding / subtracting 1 to get these
    // to match
    label sourceIndex_;

    // Component index - for Types with multiple components
    label componentIndex_;

    // Rather than store the keywords that have to be searched in a dictionary,
    // equationReader keeps its own list of keywords, and the dictLookupIndex
    // is the index in this list.  This is only applicable if the sourceType is
    // of type sldictSource
    label dictLookupIndex_;

    // The operation to be performed (+ - sin exp min, etc...)
    operationType operation_;

    // A pointer to the scalarField source data retrieval function
    mutable const scalarField&
    (Foam::equationReader::*getSourceScalarFieldFunction_)
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    // A pointer to the scalarField operation to be performed
    mutable void (Foam::equationReader::*opScalarFieldFunction_)
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storageIndex,
        scalarField& x,
        const scalarField& source
    ) const;

    // A pointer to the scalar source data retrieval function
    mutable scalar (Foam::equationReader::*getSourceScalarFunction_)
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    // A pointer to the scalar operation to be performed
    mutable void (Foam::equationReader::*opScalarFunction_)
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storageIndex,
        scalar& x,
        scalar source
    ) const;

    // A pointer to the dimensions source data retrieval function
    mutable dimensionSet (Foam::equationReader::*getSourceDimsFunction_)
    (
        const equationReader * eqnReader,
        const label equationIndex,
        const label equationOperationIndex,
        const label maxStoreIndex,
        const label storageOffset
    ) const;

    // A pointer to the dimension operation to be performed
    mutable void (Foam::equationReader::*opDimsFunction_)
    (
        const equationReader * eqnReader,
        const label index,
        const label i,
        const label storageOffset,
        label& storageIndex,
        dimensionSet& xDims,
        dimensionSet sourceDims
    ) const;


public:

    // Static data members

        static const char* const typeName;

    // Constructors

        //- Construct null
        equationOperation();

        //- Construct copy
        equationOperation(const equationOperation& eqop);

        //- Construct from components
        equationOperation
        (
            sourceTypeEnum source,
            label sourceIndex,
            label componentIndex,
            label dictLookupIndex,
            operationType operation,
            const scalarField&
            (Foam::equationReader::*getSourceScalarFieldFunction_)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                const label
            ) const = NULL,
            void (Foam::equationReader::*opScalarFieldFunction_)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                label&,
                scalarField&,
                const scalarField&
            ) const = NULL,
            scalar (Foam::equationReader::*getSourceScalarFunction)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                const label
            ) const = NULL,
            void (Foam::equationReader::*opScalarFunction)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                label&,
                scalar&,
                scalar
            ) const = NULL,
            dimensionSet (Foam::equationReader::*getSourceDimsFunction)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                const label
            ) const = NULL,
            void (Foam::equationReader::*opDimsFunction)
            (
                const equationReader *,
                const label,
                const label,
                const label,
                label&,
                dimensionSet&,
                dimensionSet
            ) const = NULL
        );


    // Destructor
    ~equationOperation();

    // Member functions

        // Access

            //- Const access to source
            inline const sourceTypeEnum& sourceType() const;

            //- Access to source
            inline sourceTypeEnum& sourceType();

            //- Const access to source index
            inline const label& sourceIndex() const;

            //- Access to source index
            inline label& sourceIndex();

            //- Const access to componentIndex
            inline const label& componentIndex() const;

            //- Access to componentIndex
            inline label& componentIndex();

            //- Const access to dictionary lookup name index
            inline const label& dictLookupIndex() const;

            //- Access to dictionary lookup name index
            inline label& dictLookupIndex();

            //- Const access to operation
            inline const operationType& operation() const;

            //- Access to operation
            inline operationType& operation();

        // Function pointers

            //- Assign the source scalarField function
            void assignSourceScalarFieldFunction
            (
                const scalarField&
                (Foam::equationReader::*getSourceScalarFieldFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    const label
                ) const
            ) const;

            //- Assign the operation scalarField function
            void assignOpScalarFieldFunction
            (
                void (Foam::equationReader::*opScalarFieldFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    label&,
                    scalarField&,
                    const scalarField&
                ) const
            ) const;

            //- Assign the source scalar function
            void assignSourceScalarFunction
            (
                scalar (Foam::equationReader::*getSourceScalarFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    const label
                ) const
            ) const;

            //- Assign the operation scalar function
            void assignOpScalarFunction
            (
                void (Foam::equationReader::*opScalarFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    label&,
                    scalar&,
                    scalar
                ) const
            ) const;

            //- Assign the source dimensions function
            void assignSourceDimsFunction
            (
                dimensionSet (Foam::equationReader::*getSourceDimsFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    const label
                ) const
            ) const;

            //- Assign the operation dimensions function
            void assignOpDimsFunction
            (
                void (Foam::equationReader::*opDimsFunction)
                (
                    const equationReader *,
                    const label,
                    const label,
                    const label,
                    label&,
                    dimensionSet&,
                    dimensionSet
                ) const
            ) const;

            //- Call the getSourceScalarField function
            const scalarField& getSourceScalarFieldFunction
            (
                const equationReader * eqnReader,
                const label equationIndex,
                const label equationOperationIndex,
                const label maxStoreIndex,
                const label storageOffset
            ) const;

            //- Call the operation scalarField function
            void opScalarFieldFunction
            (
                const equationReader * eqnReader,
                const label index,
                const label i,
                const label storageOffset,
                label& storageIndex,
                scalarField& x,
                const scalarField& source
            ) const;

            //- Call the getSourceScalar function
            scalar getSourceScalarFunction
            (
                const equationReader * eqnReader,
                const label equationIndex,
                const label equationOperationIndex,
                const label maxStoreIndex,
                const label storageOffset
            ) const;

            //- Call the operation scalar function
            void opScalarFunction
            (
                const equationReader * eqnReader,
                const label index,
                const label i,
                const label storageOffset,
                label& storageIndex,
                scalar& x,
                scalar source
            ) const;

            //- Call the getSourceDims function
            dimensionSet getSourceDimsFunction
            (
                const equationReader * eqnReader,
                const label equationIndex,
                const label equationOperationIndex,
                const label maxStoreIndex,
                const label storageOffset
            ) const;

            void opDimsFunction
            (
                const equationReader * eqnReader,
                const label index,
                const label i,
                const label storageOffset,
                label& storageIndex,
                dimensionSet& xDims,
                dimensionSet sourceDims
            ) const;

        // Convenience

            //- Look up operation number, given a word
            static operationType findOp(const word& opName);

            //- Look up operation name
            static word opName(const operationType& op);

            //- Look up sourceType name
            static word sourceName(const sourceTypeEnum& st);

    // Operators

        void operator=(equationOperation&);


    // Friend IOstream Operators

        friend Istream& operator>>(Istream&, equationOperation&);
        friend Ostream& operator<<(Ostream&, const equationOperation&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationOperationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
