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
    Foam::pisoFaControl

Description
    Specialization of the pimpleFaControl class for PISO control.

\*---------------------------------------------------------------------------*/

#ifndef pisoFaControl_H
#define pisoFaControl_H

#include "pimpleFaControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pisoFaControl Declaration
\*---------------------------------------------------------------------------*/

class pisoFaControl
:
    public pimpleFaControl
{
    // Private member functions

        //- Disallow default bitwise copy construct
        pisoFaControl(const pisoFaControl&);

        //- Disallow default bitwise assignment
        void operator=(const pisoFaControl&);


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pisoFaControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        pisoFaControl(faMesh& mesh, const word& dictName="PISO");


    //- Destructor
    virtual ~pisoFaControl();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
