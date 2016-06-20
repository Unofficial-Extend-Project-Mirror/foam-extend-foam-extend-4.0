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

#ifndef makeCoalParcelSurfaceReactionModels_H
#define makeCoalParcelSurfaceReactionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "ReactingMultiphaseCloudTemplate.H"

#include "NoSurfaceReaction.H"
#include "COxidationDiffusionLimitedRate.H"
#include "COxidationKineticDiffusionLimitedRate.H"
#include "COxidationMurphyShaddix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeCoalSurfaceReactionModels(ParcelType)                             \
                                                                              \
    makeCoalSurfaceReactionModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
    makeCoalSurfaceReactionModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
    makeCoalSurfaceReactionModelThermoType                                    \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeCoalSurfaceReactionModelThermoType(ParcelType, ThermoType)        \
                                                                              \
    makeSurfaceReactionModel                                                  \
    (                                                                         \
        ReactingMultiphaseCloud<ParcelType<ThermoType> >                      \
    );                                                                        \
                                                                              \
    makeSurfaceReactionModelThermoType                                        \
    (                                                                         \
        NoSurfaceReaction,                                                    \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoType                                        \
    (                                                                         \
        COxidationDiffusionLimitedRate,                                       \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoType                                        \
    (                                                                         \
        COxidationKineticDiffusionLimitedRate,                                \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeSurfaceReactionModelThermoType                                        \
    (                                                                         \
        COxidationMurphyShaddix,                                              \
        ReactingMultiphaseCloud,                                              \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
