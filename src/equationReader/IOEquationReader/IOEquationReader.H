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
    Foam::IOEquationReader

Description
    IOEquationReader is an equationReader with IOobject functionality for easy
    input and output.  Originally created so equations can be written out at
    every timestep if desired.

SourceFiles
    IOEquationReader.C
    IOEquationReaderIO.C

\*---------------------------------------------------------------------------*/

#ifndef IOEquationReader_H
#define IOEquationReader_H

#include "equationReader.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class IOEquationReader Declaration
\*---------------------------------------------------------------------------*/

class IOEquationReader
:
    public regIOobject,
    public equationReader
{

    // Private member data

        //- Output flag - true = include data source info in output file
        bool showDataSourceInfo_;

public:

    TypeName("equationReader");


    // Constructors

        //- Construct given an IOobject.  When true, showDataSourceInfo causes
        //  IOequationReader to include the dataSource info in its output file,
        //  even though it cannot read this data.
        IOEquationReader
        (
            const IOobject&,
            const bool showDataSourceInfo = false
        );


    // Destructor

        virtual ~IOEquationReader();


    // Member functions

        //- Access showDataSourceInfo flag
        inline bool& showDataSourceInfo()
        {
            return showDataSourceInfo_;
        }

        //- Name function is needed to disambiguate those inherited
        //  from regIOobject and dictionary
        const word& name() const;

        //- ReadData function required for regIOobject read operation
        virtual bool readData(Istream&);

        //- WriteData function required for regIOobject write operation
        virtual bool writeData(Ostream&) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
