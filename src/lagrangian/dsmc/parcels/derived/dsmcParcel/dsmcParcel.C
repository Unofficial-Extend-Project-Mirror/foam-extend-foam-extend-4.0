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

\*---------------------------------------------------------------------------*/

#include "dsmcParcel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(dsmcParcel, 0);
    defineParticleTypeNameAndDebug(dsmcParcel, 0);
    defineParcelTypeNameAndDebug(dsmcParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dsmcParcel::dsmcParcel
(
    DsmcCloud<dsmcParcel>& owner,
    const vector& position,
    const vector& U,
    const scalar Ei,
    const label celli,
    const label typeId
)
:
    DsmcParcel<dsmcParcel>
    (
        owner,
        position,
        U,
        Ei,
        celli,
        typeId
    )
{}


Foam::dsmcParcel::dsmcParcel
(
    const Cloud<dsmcParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    DsmcParcel<dsmcParcel>(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::dsmcParcel::~dsmcParcel()
{}


// ************************************************************************* //
