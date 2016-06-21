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
    Foam::pisoControl

Description
    Specialization of the pimpleControl class for PISO control.

\*---------------------------------------------------------------------------*/

#ifndef pisoControl_H
#define pisoControl_H

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class pisoControl Declaration
\*---------------------------------------------------------------------------*/

class pisoControl
:
    public pimpleControl
{
    // Private member functions

        //- Disallow default bitwise copy construct
        pisoControl(const pisoControl&);

        //- Disallow default bitwise assignment
        void operator=(const pisoControl&);


public:

    // Static Data Members

        //- Run-time type information
        TypeName("pisoControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        pisoControl(fvMesh& mesh, const word& dictName="PISO");


    //- Destructor
    virtual ~pisoControl();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
