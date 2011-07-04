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

#include "equationOperationList.H"
//#include "dimensionedScalar.H"
#include "scalarList.H"
#include "dimensionSet.H"

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
:
    public equationOperationList
{

    // Private data

        word equationName_;
        
        //- raw text read from the dictionary
        string rawText_;

        //- Result of most recent evaluate() or update()
        dimensionedScalar lastResult_;

        //- Constants appearing in the equation are stored here
//        scalarList * internalScalars_;
        
        //- Names of variables that need to be read from a dictionary are
        //  stored here
//        wordList dictLookups_;
        
        //- Override dimensions - if read from a dictionary with dimensions,
        //  eg:  nu      nu [0 2 -1 0 0 0 0] "sin(theta)";
        //  the dimensionedScalar resulting from evaluate() is given these
        //  dimensions
        dimensionSet overrideDimensions_;
        
        //- true if there is a dimension override
        bool changeDimensions_;

public:

    // Static data members
    
        static const char* const typeName;


    // Constructors

        //- Construct null
        equation();
        
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
            inline const word& equationName() const;

            //- Change the equation name
            inline word& equationName();
            
            //- Equation text
            inline const string& rawText() const;
            
            //- Set the equation text
            inline string& rawText();

            //- Last result
            inline const dimensionedScalar& lastResult() const;
            
            //- Set the last result
            inline dimensionedScalar& lastResult();

            //- Easy access to operations
            inline const equationOperationList& ops() const;

            inline equationOperationList& ops();
            
            //- Locally stored constants
//            inline const scalarList * internalScalars() const;

            //- Set the locally stored constants
//            inline scalarList * internalScalars();
            
            //- Dictionary variable lookup names
//            inline const wordList& dictLookups() const;

            //- Set the dictionary variable lookup names
//            inline wordList& dictLookups();
            
            //- Dimension override
            inline const dimensionSet& overrideDimensions() const;

            //- Set the dimension override
            inline dimensionSet& overrideDimensions();

            //- changeDimensions flag
            inline const bool& changeDimensions() const;

            //- Set the changeDimensions flag
            inline bool& changeDimensions();

    // Operators
        // Copy only the header info - rawText, equationName, dimensionOverride
        // and changeDimensions
        void operator=(equation&);
      
    // Friend IOstream Operators

        friend Istream& operator>>(Istream&, equation&);
        friend Ostream& operator<<(Ostream&, const equation&);
      
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "equationI.H"

#endif
