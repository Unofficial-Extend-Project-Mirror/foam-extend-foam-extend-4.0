/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     3.2
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

#include "eConstThermo.H"

#include "hConstThermo.H"
#include "janafThermo.H"
#include "specieThermo.H"

#include "constTransport.H"
#include "sutherlandTransport.H"

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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
