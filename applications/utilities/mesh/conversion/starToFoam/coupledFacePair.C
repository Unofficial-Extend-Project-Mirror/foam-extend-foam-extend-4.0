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

Description
    Data associated with a pair of coupled faces.
    1 represents integral match; all other number are arbitrary matches

\*---------------------------------------------------------------------------*/

#include "coupledFacePair.H"
#include "error.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
coupledFacePair::coupledFacePair
(
    const label coupleNo,
    const label mC, const label mF,
    const label sC, const label sF,
    const label integral
)
:
    coupleID_(coupleNo),
    masterCellID_(mC),
    masterFaceID_(mF),
    slaveCellID_(sC),
    slaveFaceID_(sF),
    integralMatch_(integral == 1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Ostream& operator<<(Ostream& os, const coupledFacePair& c)
{
    os  << "Master cell: " << c.masterCellID_
        << " face: " << c.masterFaceID_ << endl
        << "Slave cell: " << c.slaveCellID_
        << " face: " << c.slaveFaceID_ << endl
        << "Integral: " << c.integralMatch_ << endl;

    return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
