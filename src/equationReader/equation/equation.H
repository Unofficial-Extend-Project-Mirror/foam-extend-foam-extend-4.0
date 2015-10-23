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
    Foam::equation

Description
    An equation, read from a dictionary.  Effectively just a container class
    holding all data associated with an individual equation.  Functions are
    implemented in masterDictionary.

SourceFiles
    equationI.H
    equation.C
    equationIO.C

Author
    David L. F. Gaden

\*---------------------------------------------------------------------------*/

#ifndef equation_H
#define equation_H

//#include "equationOperationList.H"
//#include "dimensionedScalar.H"
#include "scalarField.H"
#include "dimensionSet.H"
#include "equationOperation.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

class equation;

// Friend IOstream Operators

Istream& operator>>(Istream&, equation&);
Ostream& operator<<(Ostream&, const equation&);

/*---------------------------------------------------------------------------*\
                        Class equation Declaration
\*---------------------------------------------------------------------------*/

class equation
{

    // Private data

        word equationName_;

        //- Operation list
        mutable PtrList<equationOperation> ops_;

        //- raw text read from the dictionary
        mutable string rawText_;

        //- Result of most recent evaluate() or update()
        mutable dimensionedScalar lastResult_;

        //- Override dimensions - if read from a dictionary with dimensions,
        //  eg:  nu      nu [0 2 -1 0 0 0 0] "sin(theta)";
        //  the dimensionedScalar resulting from evaluate() is given these
        //  dimensions
        mutable dimensionSet overrideDimensions_;

        //- true if there is a dimension override
        mutable bool changeDimensions_;

        //- Maximum field size of output, indexed by geoIndex
        mutable labelList maxFieldSizes_;

    // Private member functions

        //- Disallow default bitwise assignment
        void operator=(const equation&);

public:

    // Static data members

        static const char* const typeName;


    // Constructors

        //- Construct null
        equation();

        //- Construct copy
        equation(const equation&);

        //- Construct from Istream with optional name
        equation(Istream& is, const word& name = word::null);

        //- Construct from components
        equation
        (
            const word& equationName,
            const string& rawText,
            const dimensionSet& overrideDimensions = dimless,
            const bool& changeDimensions = false
        );

    // Destructor
    ~equation();

    // Member functions

        // Access

            //- Equation name
            inline const word& name() const;

            //- Size of operation list
            inline label size() const;

            //- setSize of operation list
            inline void setSize(const label newSize) const;

            //- set an element of operation list
            inline autoPtr<equationOperation> set
            (
                const label elementIndex,
                equationOperation * newOperation
            ) const;

            //- Change the equation name
            inline word& name();

            //- Equation text
            inline const string& rawText() const;

            //- Set equation text
            inline void setRawText(const string& newRawText) const;

            //- Last result
            inline const dimensionedScalar& lastResult() const;

            //- Set the last result
            inline void setLastResult
            (
                const dimensionedScalar& newResult
            ) const;

            //- Set the value of the last result
            inline void setLastResult(const word& newName) const;

            //- Set the dimensions of the last result
            inline void setLastResult(const dimensionSet& newDims) const;

            //- Set the value of the last result
            inline void setLastResult(const scalar& newScalar) const;

            //- Dimension override
            inline const dimensionSet& overrideDimensions() const;

            //- Set override dimensions
            inline void setOverrideDimensions
            (
                const dimensionSet& newDims
            ) const;

            //- changeDimensions flag
            inline const bool& changeDimensions() const;

            //- Set the changeDimensions flag
            inline void setChangeDimensions(bool newFlag) const;

            //- Return maxFieldSizes
            inline const labelList& maxFieldSizes() const;

            //- Set maxFieldSizes
            inline void setMaxFieldSizes
            (
                const labelList& newSizes
            ) const;

    // Delete the operation list
    void clear() const;

    // Operators

        //- Copy only the header info - rawText, equationName,
        // dimensionOverride, and changeDimensions
        void operator=(equation&);

        //- Access to ops_
        inline equationOperation& operator[](const label) const;

    // Friend IOstream Operators

        friend Istream& operator>>(Istream&, equation&);
        friend Ostream& operator<<(Ostream&, const equation&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
