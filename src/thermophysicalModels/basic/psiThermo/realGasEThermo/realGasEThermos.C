/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | 
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
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    sutherlandTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    pengRobinson
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    aungierRedlichKwong
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    redlichKwong
);

makeBasicRealGasThermo
(
    realGasEThermo,
    pureMixture,
    constTransport,
    realGasSpecieThermo,
    nasaHeatCapacityPolynomial,
    soaveRedlichKwong
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
