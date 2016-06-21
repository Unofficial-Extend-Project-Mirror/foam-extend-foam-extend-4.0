/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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
    Foam::OFstream

Description
    Output to file stream.

SourceFiles
    OFstream.C

\*---------------------------------------------------------------------------*/

#ifndef OFstream_H
#define OFstream_H

#include "OSstream.H"
#include "fileName.H"
#include "className.H"

#include <fstream>
using std::ofstream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class OFstream;

/*---------------------------------------------------------------------------*\
                       Class OFstreamAllocator Declaration
\*---------------------------------------------------------------------------*/

//- A std::ostream with ability to handle compressed files
class OFstreamAllocator
{
    friend class OFstream;

    ostream* ofPtr_;

    // Constructors

        //- Construct from pathname
        OFstreamAllocator
        (
            const fileName& pathname,
            ios_base::openmode mode,
            IOstream::compressionType compression = IOstream::UNCOMPRESSED

        );

    // Destructor

        ~OFstreamAllocator();


public:

    // Member functions

        //- Access to underlying std::ostream
        ostream& stdStream();
};


/*---------------------------------------------------------------------------*\
                           Class OFstream Declaration
\*---------------------------------------------------------------------------*/

class OFstream
:
    private OFstreamAllocator,
    public OSstream
{
    // Private data

        fileName pathname_;


public:

    // Declare name of the class and its debug switch
    ClassName("OFstream");


    // Constructors

        //- Construct from pathname
        OFstream
        (
            const fileName& pathname,
            ios_base::openmode mode = ios_base::out|ios_base::trunc,
            streamFormat format = ASCII,
            versionNumber version = currentVersion,
            compressionType compression = UNCOMPRESSED
        );


    // Destructor

        ~OFstream();


    // Member functions

        // Access

            //- Return the name of the stream
            const fileName& name() const
            {
                return pathname_;
            }

            //- Return non-const access to the name of the stream
            fileName& name()
            {
                return pathname_;
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global predefined null output stream

extern OFstream Snull;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
