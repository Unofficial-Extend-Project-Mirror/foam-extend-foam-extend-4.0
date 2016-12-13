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

Description
    Mixture instantiation

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "basicMixture.H"
#include "makeBasicMixture.H"

#include "perfectGas.H"
#include "redlichKwong.H"
#include "pengRobinson.H"
#include "aungierRedlichKwong.H"
#include "soaveRedlichKwong.H"

#include "eConstThermo.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "nasaHeatCapacityPolynomial.H"
#include "constantHeatCapacity.H"
#include "specieThermo.H"
#include "realGasSpecieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "pureMixture.H"

#include "thermoPhysicsTypes.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * Private Static Data * * * * * * * * * * * * * */

makeBasicMixture
(
    pureMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    hConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    constTransport,
    eConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    eConstThermo,
    perfectGas
);

makeBasicMixture
(
    pureMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeBasicMixturePhys
(
    pureMixture,
    icoPoly8ThermoPhysics
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
