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
#include "constRealGasTransport.H"
#include "sutherlandTransport.H"
#include "sutherlandRealGasTransport.H"
#include "constantHeatCapacity.H"
#include "pureMixture.H"
#include "realGasEThermo.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/* * * * * * * * * * * * * * * private static data * * * * * * * * * * * * * */

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    pengRobinson
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    aungierRedlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    redlichKwong
);


makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constRealGasTransport,
    realGasSpecieThermo,
    constantHeatCapacity,
    soaveRedlichKwong
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
