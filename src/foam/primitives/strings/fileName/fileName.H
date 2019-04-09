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
    Foam::fileName

Description
    A class for handling file names.

    A fileName is a string of characters without whitespace or quotes.
    A fileName can be
      - constructed from a char*, a string or a word
      - concatenated by adding a '/' separator
      - decomposed into the path, name or component list
      - interrogated for type and access mode

    The string::expand() method expands environment variables, etc,

SourceFiles
    fileName.C
    fileNameIO.C

\*---------------------------------------------------------------------------*/

#ifndef fileName_H
#define fileName_H

#include "word.H"
#include "debugSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class T> class List;
typedef List<word> wordList;

// Forward declaration of friend functions and operators

class fileName;

Istream& operator>>(Istream&, fileName&);
Ostream& operator<<(Ostream&, const fileName&);


/*---------------------------------------------------------------------------*\
                          Class fileName Declaration
\*---------------------------------------------------------------------------*/

class fileName
:
    public string
{
    // Private member functions

        //- Strip invalid characters
        inline void stripInvalid();


public:

    //- Enumerations to handle file types and modes.
    enum Type
    {
        UNDEFINED,
        FILE,
        DIRECTORY,
        LINK
    };


    // Static data members

        static const char* const typeName;
        static debug::debugSwitch debug;
        static const fileName null;


    // Constructors

        //- Construct null
        inline fileName();

        //- Construct as copy
        inline fileName(const fileName&);

        //- Construct as copy of word
        inline fileName(const word&);

        //- Construct as copy of string
        inline fileName(const string&);

        //- Construct as copy of std::string
        inline fileName(const std::string&);

        //- Construct as copy of character array
        inline fileName(const char*);

        //- Construct by concatenating elements of wordList separated by '/'
        explicit fileName(const wordList&);

        //- Construct from Istream
        fileName(Istream&);


    // Member functions

        //- Is this character valid for a fileName?
        inline static bool valid(char);

        //- Cleanup file name
        //
        // * Removes repeated slashes
        //       /abc////def        -->   /abc/def
        //
        // * Removes '/./'
        //       /abc/def/./ghi/.   -->   /abc/def/./ghi
        //       abc/def/./         -->   abc/def
        //
        // * Removes '/../'
        //       /abc/def/../ghi/jkl/nmo/..   -->   /abc/ghi/jkl
        //       abc/../def/ghi/../jkl        -->   abc/../def/jkl
        //
        // * Removes trailing '/'
        //
        bool clean();

        //- Cleanup file name
        //  eg, remove repeated slashes, etc.
        fileName clean() const;


        // Interrogation

            //- Return the file type: FILE, DIRECTORY or UNDEFINED
            Type type() const;

            //- Return true if file name is absolute
            bool isAbsolute() const;

            //- Convert from relative to absolute
            fileName& toAbsolute();


        // Decomposition

            //- Return file name (part beyond last /)
            //
            //  Behaviour compared to /usr/bin/basename:
            //    Input           name()          basename
            //    -----           ------          --------
            //    "foo"           "foo"           "foo"
            //    "/foo"          "foo"           "foo"
            //    "foo/bar"       "bar"           "bar"
            //    "/foo/bar"      "bar"           "bar"
            //    "/foo/bar/"     ""              "bar"
            //
            word name() const;

            //- Return file name (part beyond last /), subsitute for FOAM_CASE
            string caseName() const;

            //- Return file name, optionally without extension
            word name(const bool noExt) const;

            //- Return directory path name (part before last /)
            //
            //  Behaviour compared to /usr/bin/dirname:
            //    input           path()          dirname
            //    -----           ------          -------
            //    "foo"           "."             "."
            //    "/foo"          "/"             "foo"
            //    "foo/bar"       "foo"           "foo"
            //    "/foo/bar"      "/foo"          "/foo"
            //    "/foo/bar/"     "/foo/bar/"     "/foo"
            //
            fileName path() const;

            //- Return file name without extension (part before last .)
            fileName lessExt() const;

            //- Return file name extension (part after last .)
            word ext() const;

            //- Return path components as wordList
            //
            //  Behaviour:
            //    Input           components()
            //    -----           ------
            //    "foo"           1("foo")
            //    "/foo"          1("foo")
            //    "foo/bar"       2("foo", "bar")
            //    "/foo/bar"      2("foo", "bar")
            //    "/foo/bar/"     2("foo", "bar")
            wordList components(const char delimiter='/') const;

            //- Return a single component of the path
            word component(const size_type, const char delimiter='/') const;


    // Member operators

        // Assignment

            void operator=(const fileName&);
            void operator=(const word&);
            void operator=(const string&);
            void operator=(const std::string&);
            void operator=(const char*);


    // IOstream operators

        friend Istream& operator>>(Istream&, fileName&);
        friend Ostream& operator<<(Ostream&, const fileName&);
};


//- Assemble words and fileNames as pathnames by adding a '/' separator
fileName operator/(const string&, const string&);


//- Recursively search the given directory for the file
//  returning the path relative to the directory or
//  fileName::null if not found
fileName search(const word& file, const fileName& directory);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "fileNameI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
