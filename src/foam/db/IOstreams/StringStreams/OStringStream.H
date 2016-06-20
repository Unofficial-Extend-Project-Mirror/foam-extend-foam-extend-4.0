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
    Foam::OStringStream

Description
    Output to memory buffer stream.

SourceFiles
    StringStreamsPrint.C

\*---------------------------------------------------------------------------*/

#ifndef OStringStream_H
#define OStringStream_H

#include "OSstream.H"
#include <sstream>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class OStringStream Declaration
\*---------------------------------------------------------------------------*/

class OStringStream
:
    public OSstream
{

public:

    // Constructors

        //- Construct and set stream status
        OStringStream
        (
            streamFormat format = ASCII,
            versionNumber version = currentVersion
        )
        :
            OSstream
            (
               *(new std::ostringstream()),
                "OStringStream.sinkFile",
                format,
                version
            )
        {}

        //- Construct as copy
        OStringStream(const OStringStream& oss)
        :
            OSstream
            (
               *(
                    new std::ostringstream
                    (
                        dynamic_cast<const std::ostringstream&>
                        (
                            oss.stream()
                        ).str()
                    )
                ),
                oss.name(),
                oss.format(),
                oss.version()
            )
        {}


    // Destructor

        ~OStringStream()
        {
            delete &dynamic_cast<std::ostringstream&>(stream());
        }


    // Member functions

        // Access

            //- Return the string
            string str() const
            {
                return dynamic_cast<const std::ostringstream&>(stream()).str();
            }


        // Edit

            //- Clear the OStringStream
            void rewind()
            {
#               if __GNUC__ < 4 && __GNUC_MINOR__ < 4
                stream().rdbuf()->pubsetbuf(" ", 1);
#               else
                stream().rdbuf()->pubseekpos(0);
#               endif
            }


        // Print

            //- Print description to Ostream
            void print(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
