/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
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
    Foam::dynamicCodeContext

Description
    Encapsulation of dynamic code dictionaries

SourceFiles
    dynamicCodeContext.C

\*---------------------------------------------------------------------------*/

#ifndef dynamicCodeContext_H
#define dynamicCodeContext_H

#include "dictionary.H"
#include "SHA1Digest.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class dynamicCodeContext Declaration
\*---------------------------------------------------------------------------*/

class dynamicCodeContext
{
    // Private data

        //- The parent dictionary context
        const dictionary& dict_;

        //- Mandatory "code" entry
        string code_;

        //- Optional "localCode" entry
        string localCode_;

        //- Optional "codeInclude" entry
        string include_;

        //- Optional "codeOptions" entry
        string options_;

        //- Optional "codeLib" entry
        string libs_;

        //- Calculated SHA1Digest
        SHA1Digest sha1_;


public:

    // Constructors

        //- Construct from a dictionary
        dynamicCodeContext(const dictionary&);


    // Member functions

        //- Return the parent dictionary context
        const dictionary& dict() const
        {
            return dict_;
        }

        //- Return the code-includes
        const string& include() const
        {
            return include_;
        }

        //- Return the code-options
        const string& options() const
        {
            return options_;
        }

        //- Return the code-libs
        const string& libs() const
        {
            return libs_;
        }

        //- Return the code
        const string& code() const
        {
            return code_;
        }

        //- Return the local (file-scope) code
        const string& localCode() const
        {
            return localCode_;
        }

        //- Return SHA1 digest calculated from include, options, code
        const SHA1Digest& sha1() const
        {
            return sha1_;
        }

        //- Helper: add \#line directive
        static void addLineDirective
        (
            string&,
            const label lineNum,
            const fileName& name
        );

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
