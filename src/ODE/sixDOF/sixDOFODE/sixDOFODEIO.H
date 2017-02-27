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
    sixDOFODEIO

Description
    Class used to control input/output for sixDOFODEIO. Need separate classes in
    order to be able to have run-time selection and automatic I/O.

Author
    Vuko Vukcevic, FSB Zagreb.  All rights reserved.

SourceFiles
    sixDOFODEIO.C

\*---------------------------------------------------------------------------*/

#ifndef sixDOFODEIO_H
#define sixDOFODEIO_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sixDOFODE;

/*---------------------------------------------------------------------------*\
                         Class sixDOFODEIO Declaration
\*---------------------------------------------------------------------------*/

class sixDOFODEIO
:
    public IOdictionary
{
    // Private data

        //- Const reference to the underlying sixDOFODE object
        const sixDOFODE& sixDOF_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        sixDOFODEIO(const sixDOFODEIO&);

        //- Disallow default bitwise assignment
        void operator=(const sixDOFODEIO&);


public:

    // Constructors

        //- Construct from dictionary and sixDOFODE object
        sixDOFODEIO(const IOobject& io, const sixDOFODE& sixDOF);


    // Destructor

        virtual ~sixDOFODEIO();


        // Write control

            //- writeData member function controlling output
            virtual bool writeData(Ostream& os) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
