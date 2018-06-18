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

Author
Christian Lucas
Institut für Thermodynamik
Technische Universität Braunschweig
Germany

\*---------------------------------------------------------------------------*/

#include "makeBasicPsiThermo.H"

#include "redlichKwong.H"
#include "pengRobinson.H"
#include "aungierRedlichKwong.H"
#include "soaveRedlichKwong.H"
#include "nasaHeatCapacityPolynomial.H"
#include "realGasSpecieThermo.H"
#include "constTransport.H"
#include "sutherlandTransport.H"
#include "constantHeatCapacity.H"
#include "pureMixture.H"
#include "realGasHThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);


makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasHThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
