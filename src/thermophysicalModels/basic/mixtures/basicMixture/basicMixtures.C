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

Description
    Mixture instantiation

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "basicMixture.H"
#include "makeBasicMixture.H"

#include "perfectGas.H"

#include "eConstThermo.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "pureMixture.H"

#include "thermoPhysicsTypes.H"

#include "addToRunTimeSelectionTable.H"

#include "redlichKwong.H"
#include "pengRobinson.H"
#include "aungierRedlichKwong.H"
#include "soaveRedlichKwong.H"
#include "realGasSpecieThermo.H"
#include "nasaHeatCapacityPolynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

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


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
