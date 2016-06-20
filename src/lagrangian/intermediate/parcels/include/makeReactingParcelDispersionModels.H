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

#ifndef makeReactingParcelDispersionModels_H
#define makeReactingParcelDispersionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "KinematicCloudTemplate.H"

#include "NoDispersion.H"
#include "GradientDispersionRAS.H"
#include "StochasticDispersionRAS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingDispersionModels(ParcelType)                              \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingDispersionModelThermoType                                     \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingDispersionModelThermoType(ParcelType, ThermoType)         \
                                                                              \
    makeDispersionModel(KinematicCloud<ParcelType<ThermoType> >);             \
                                                                              \
    defineNamedTemplateTypeNameAndDebug                                       \
    (                                                                         \
        DispersionRASModel<KinematicCloud<ParcelType<ThermoType> > >,         \
        0                                                                     \
    );                                                                        \
                                                                              \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        NoDispersion,                                                         \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        GradientDispersionRAS,                                                \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeDispersionModelThermoType                                             \
    (                                                                         \
        StochasticDispersionRAS,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
