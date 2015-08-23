/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     |
    \\  /    A nd           | For copyright notice see file Copyright
     \\/     M anipulation  |
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
#include "constRealGasTransport.H"
#include "sutherlandRealGasTransport.H"

#include "pureMixture.H"

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

//CL: Real Gas Mixtures
makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealFluidMixture
(
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
