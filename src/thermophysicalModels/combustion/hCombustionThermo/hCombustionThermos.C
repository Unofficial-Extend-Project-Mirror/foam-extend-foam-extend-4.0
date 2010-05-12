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

#include "hCombustionThermo.H"
#include "hMixtureThermo.H"

#include "makeCombustionThermo.H"
#include "addToRunTimeSelectionTable.H"

#include "perfectGas.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

#include "homogeneousMixture.H"
#include "inhomogeneousMixture.H"
#include "veryInhomogeneousMixture.H"
#include "dieselMixture.H"
#include "multiComponentMixture.H"
#include "reactingMixture.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    homogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    inhomogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    veryInhomogeneousMixture,
    constTransport,
    hConstThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    homogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    inhomogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    veryInhomogeneousMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);

makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    dieselMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);
   
makeCombustionThermo
(
    hCombustionThermo,
    hMixtureThermo,
    multiComponentMixture,
    sutherlandTransport,
    janafThermo,
    perfectGas
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Reaction thermo

defineTemplateTypeNameAndDebug(hMixtureThermo<reactingMixture>, 0)

typedef hMixtureThermo<reactingMixture> hMixtureThermoReactingMixture;

addToRunTimeSelectionTable
(
    hCombustionThermo,
    hMixtureThermoReactingMixture,
    fvMesh
);

addToRunTimeSelectionTable
(
    basicThermo,
    hMixtureThermoReactingMixture,
    fvMesh
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
