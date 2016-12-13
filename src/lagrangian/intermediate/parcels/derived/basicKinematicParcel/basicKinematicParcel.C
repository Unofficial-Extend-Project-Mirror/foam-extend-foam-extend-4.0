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

#include "basicKinematicParcel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicKinematicParcel, 0);
    defineParticleTypeNameAndDebug(basicKinematicParcel, 0);
    defineParcelTypeNameAndDebug(basicKinematicParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicKinematicParcel::basicKinematicParcel
(
    KinematicCloud<basicKinematicParcel>& owner,
    const vector& position,
    const label cellI
)
:
    KinematicParcel<basicKinematicParcel>(owner, position, cellI)
{}


Foam::basicKinematicParcel::basicKinematicParcel
(
    KinematicCloud<basicKinematicParcel>& owner,
    const vector& position,
    const label cellI,
    const label typeId,
    const scalar nParticle0,
    const scalar d0,
    const vector& U0,
    const constantProperties& constProps
)
:
    KinematicParcel<basicKinematicParcel>
    (
        owner,
        position,
        cellI,
        typeId,
        nParticle0,
        d0,
        U0,
        constProps
    )
{}


Foam::basicKinematicParcel::basicKinematicParcel
(
    const Cloud<basicKinematicParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    KinematicParcel<basicKinematicParcel>(cloud, is, readFields)
{}


Foam::basicKinematicParcel::basicKinematicParcel
(
    const basicKinematicParcel& p
)
:
    KinematicParcel<basicKinematicParcel>(p)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::basicKinematicParcel::~basicKinematicParcel()
{}


// ************************************************************************* //
