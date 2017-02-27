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

\*---------------------------------------------------------------------------*/

#include "sixDOFODEIO.H"
#include "sixDOFODE.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDOFODEIO::sixDOFODEIO(const IOobject& io, const sixDOFODE& sixDOF)
:
    IOdictionary(io),
    sixDOF_(sixDOF)
{
    // Note: parameter sixDOF is incomplete here, must not call its member
    // functions here.
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::sixDOFODEIO::~sixDOFODEIO()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::sixDOFODEIO::writeData(Ostream& os) const
{
    return sixDOF_.writeData(os);
}


// ************************************************************************* //
