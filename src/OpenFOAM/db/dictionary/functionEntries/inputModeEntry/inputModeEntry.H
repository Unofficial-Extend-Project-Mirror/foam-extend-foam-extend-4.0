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
    Foam::functionEntries::inputModeEntry

Description
    Specify the input mode when reading dictionaries, expects
    a single word to follow.

    An example of @c \#inputMode directive:
    @verbatim
        #inputMode merge
    @endverbatim

    The possible input modes:
    @param merge      merge sub-dictionaries when possible
    @param overwrite  keep last entry and silently remove previous ones
    @param error      flag duplicate entry as an error
    @param default    currently the same as error

SourceFiles
    inputModeEntry.C

\*---------------------------------------------------------------------------*/

#ifndef inputModeEntry_H
#define inputModeEntry_H

#include "functionEntry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionEntries
{

/*---------------------------------------------------------------------------*\
                        Class inputModeEntry Declaration
\*---------------------------------------------------------------------------*/

class inputModeEntry
:
    public functionEntry
{
        //- input mode options
        enum inputMode
        {
            imError = 0,
            imMerge = 0x1,
            imOverwrite = 0x2
        };

        //- current input mode
        static label mode_;


    // Private Member Functions

        //- Read the mode as a word and set enum appropriately
        static void setMode(Istream&);

        //- Disallow default bitwise copy construct
        inputModeEntry(const inputModeEntry&);

        //- Disallow default bitwise assignment
        void operator=(const inputModeEntry&);


public:

    //- Runtime type information
    TypeName("inputMode");


    // Member Functions

        //- Execute the functionEntry in a sub-dict context
        static bool execute
        (
            dictionary& parentDict,
            Istream&
        );

        //- Reset the inputMode to 'default'
        static void clear();

        //- Return true if the inputMode is 'merge'
        static bool merge();

        //- Return true if the inputMode is 'overwrite'
        static bool overwrite();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionEntries
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
