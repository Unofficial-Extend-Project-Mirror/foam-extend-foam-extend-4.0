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
    Foam::primitiveEntry

Description
    A keyword and a list of tokens is a 'primitiveEntry'.
    An primitiveEntry can be read, written and printed, and the types and
    values of its tokens analysed.

    A primitiveEntry is a high-level building block for data description. It
    is a front-end for the token parser. A list of entries can be used as a
    set of keyword syntax elements, for example.

SourceFiles
    primitiveEntry.C
    primitiveEntryIO.C

\*---------------------------------------------------------------------------*/

#ifndef primitiveEntry_H
#define primitiveEntry_H

#include "IStringStream.H"
#include "OStringStream.H"

#include "entry.H"
#include "ITstream.H"
#include "InfoProxy.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;

/*---------------------------------------------------------------------------*\
                           Class primitiveEntry Declaration
\*---------------------------------------------------------------------------*/

class primitiveEntry
:
    public entry,
    public ITstream
{

public:

    // Private member functions

        //- Append the given token to this entry
        void append
        (
            const token& currToken,
            const dictionary& dict,
            Istream& is
        );

        //- Append the given tokens starting at the current tokenIndex
        void append(const tokenList& varTokens);

        //- Expand the given variable (keyword starts with $)
        bool expandVariable(const word& keyword, const dictionary& dict);

        //- Expand the given function (keyword starts with #)
        bool expandFunction
        (
            const word& keyword,
            const dictionary& dict,
            Istream& is
        );

        //- Read tokens from the given stream
        bool read(const dictionary& dict, Istream&);

        //- Read the complete entry from the given stream
        void readEntry(const dictionary& dict, Istream&);

        //- Insert the given tokens at token i
        void insert(const tokenList& varTokens, const label i);


public:

    // Constructors

        //- Construct from keyword and a Istream
        primitiveEntry(const word& keyword, Istream&);

        //- Construct from keyword, parent dictionary and a Istream
        primitiveEntry(const word& keyword, const dictionary&, Istream&);

        //- Construct from keyword and a ITstream
        primitiveEntry(const word& keyword, const ITstream&);

        //- Construct from keyword and a token
        primitiveEntry(const word&, const token&);

        //- Construct from keyword and a tokenList
        primitiveEntry(const word&, const tokenList&);

        //- Construct from keyword and a T
        template<class T>
        primitiveEntry(const word&, const T&);

        autoPtr<entry> clone(const dictionary&) const
        {
            return autoPtr<entry>(new primitiveEntry(*this));
        }


    // Member functions

        //- Return the dictionary name
        const fileName& name() const
        {
            return ITstream::name();
        }

        //- Return the dictionary name
        fileName& name()
        {
            return ITstream::name();
        }

        //- Return line number of first token in dictionary
        label startLineNumber() const;

        //- Return line number of last token in dictionary
        label endLineNumber() const;

        //- Return true because this entry is a stream
        bool isStream() const
        {
            return true;
        }

        //- Return token stream if this entry is a primitive entry
        ITstream& stream() const;

        //- This entry is not a dictionary,
        //  calling this function generates a FatalError
        const dictionary& dict() const;

        //- This entry is not a dictionary,
        //  calling this function generates a FatalError
        dictionary& dict();

        // Write
        void write(Ostream&) const;

        //- Return info proxy.
        //  Used to print token information to a stream
        InfoProxy<primitiveEntry> info() const
        {
            return *this;
        }
};


template<>
Ostream& operator<<(Ostream& os, const InfoProxy<primitiveEntry>& ip);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#   include "primitiveEntryTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
