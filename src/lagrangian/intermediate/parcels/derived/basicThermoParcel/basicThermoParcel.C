/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "basicThermoParcel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(basicThermoParcel, 0);
    defineParticleTypeNameAndDebug(basicThermoParcel, 0);
    defineParcelTypeNameAndDebug(basicThermoParcel, 0);
};


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermoParcel::basicThermoParcel
(
    ThermoCloud<basicThermoParcel>& owner,
    const label typeId,
    const vector position,
    const label celli,
    const scalar d0,
    const vector U0,
    const scalar nParticle0,
    const constantProperties& constProps
)
:
    ThermoParcel<basicThermoParcel>
    (
        owner,
        typeId,
        position,
        celli,
        d0,
        U0,
        nParticle0,
        constProps
    )
{}


Foam::basicThermoParcel::basicThermoParcel
(
    const Cloud<basicThermoParcel>& cloud,
    Istream& is,
    bool readFields
)
:
    ThermoParcel<basicThermoParcel>(cloud, is, readFields)
{}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

Foam::basicThermoParcel::~basicThermoParcel()
{}


// ************************************************************************* //
