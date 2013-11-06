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
    Foam::calcTypes::scalarMult

Description
    Multiplies base field with a scalar value.

    New field name specified by -resultName option, or automatically as:
        <baseFieldName>_scalarMult_value

    Example usage:
        scalarMult p scalarMult -value 100000 -resultName pAbs

SourceFiles
    scalarMult.C
    writeScalarMultField.C
    writeScalarMultValue.C

Author
    Dominik Christ, Wikki Ltd. All right reserved.

\*---------------------------------------------------------------------------*/

#ifndef scalarMult_H
#define scalarMult_H

#include "calcType.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace calcTypes
{

/*---------------------------------------------------------------------------*\
                          Class scalarMult Declaration
\*---------------------------------------------------------------------------*/

class scalarMult
:
    public calcType
{
    // Private data

        //- Name of base field (to scalarMult to)
        word baseFieldName_;

        //- String representation of scalar value
        string scalarMultValueStr_;

        //- Name of result field
        word resultName_;


    // Private Member Functions

        // Output

            //- Calc and output field and value scalarMult
            void writeScalarMultValues
            (
                const Time& runTime,
                const fvMesh& mesh,
                const IOobject& baseFieldHeader
            );


        //- Disallow default bitwise copy construct
        scalarMult(const scalarMult&);

        //- Disallow default bitwise assignment
        void operator=(const scalarMult&);


protected:

    // Member Functions

        // Calculation routines

            //- Initialise - typically setting static variables,
            //  e.g. command line arguments
            virtual void init();

            //- Pre-time loop calculations
            virtual void preCalc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );

            //- Time loop calculations
            virtual void calc
            (
                const argList& args,
                const Time& runTime,
                const fvMesh& mesh
            );


        // I-O

            //- Write scalarMult value
            template<class Type>
            void writeScalarMultValue
            (
                const IOobject& baseHeader,
                const string& valueStr,
                const fvMesh& mesh,
                bool& processed
            );


public:

    //- Runtime type information
    TypeName("scalarMult");


    // Constructors

        //- Construct null
        scalarMult();


    // Destructor

        virtual ~scalarMult();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace calcTypes
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "writeScalarMultValue.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
