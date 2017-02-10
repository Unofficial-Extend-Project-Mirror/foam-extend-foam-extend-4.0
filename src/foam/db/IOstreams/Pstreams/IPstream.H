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
    Foam::IPstream

Description
    Input inter-processor communications stream.

SourceFiles
    IPstream.C

\*---------------------------------------------------------------------------*/

#include "Pstream.H"

#ifndef IPstream_H
#define IPstream_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IPstream Declaration
\*---------------------------------------------------------------------------*/

class IPstream
:
    public Pstream,
    public Istream
{
    // Private data

        //- ID of sending processor
        int fromProcNo_;

        //- Message tag
        const int tag_;

        //- Communicator
        const label comm_;

        //- Message size
        label messageSize_;


    // Private member functions

        //- Check the bufferPosition_ against messageSize_ for EOF
        inline void checkEof();

        //- Read a T from the transfer buffer
        template<class T>
        inline void readFromBuffer(T&);

        //- Read data from the transfer buffer
        inline void readFromBuffer(void* data, size_t count, size_t align);


public:

    // Constructors

        //- Construct given process index to read from, optional buffer size,
        //  read format and IO version
        IPstream
        (
            const commsTypes commsType,
            const int fromProcNo,
            const label bufSize = 0,
            const int tag = Pstream::msgType(),
            const label comm = Pstream::worldComm,
            streamFormat format = BINARY,
            versionNumber version = currentVersion
        );


    //- Destructor
    ~IPstream();


    // Member functions

        // Inquiry

            //- Return flags of output stream
            ios_base::fmtflags flags() const
            {
                return ios_base::fmtflags(0);
            }


        // Read functions

            //- Read into given buffer from given processor and return the
            //  message size
            static label read
            (
                const commsTypes commsType,
                const int fromProcNo,
                char* buf,
                const std::streamsize bufSize,
                const int tag = Pstream::msgType(),
                const label communicator = 0
            );

            //- Return next token from stream
            Istream& read(token&);

            //- Read a character
            Istream& read(char&);

            //- Read a word
            Istream& read(word&);

            // Read a string (including enclosing double-quotes)
            Istream& read(string&);

            //- Read a label
            Istream& read(label&);

            //- Read a floatScalar
            Istream& read(floatScalar&);

            //- Read a doubleScalar
            Istream& read(doubleScalar&);

            //- Read a longDoubleScalar
            Istream& read(longDoubleScalar&);

            //- Read binary block
            Istream& read(char*, std::streamsize);

            //- Rewind and return the stream so that it may be read again
            Istream& rewind();


        // Edit

            //- Set flags of stream
            ios_base::fmtflags flags(const ios_base::fmtflags)
            {
                return ios_base::fmtflags(0);
            }


        // Print

            //- Print description of IOstream to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
