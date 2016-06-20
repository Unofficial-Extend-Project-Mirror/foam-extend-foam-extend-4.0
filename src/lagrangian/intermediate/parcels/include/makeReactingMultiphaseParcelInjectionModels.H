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

#ifndef makeReactingMultiphaseParcelInjectionModels_H
#define makeReactingMultiphaseParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "thermoPhysicsTypes.H"
#include "KinematicCloudTemplate.H"

#include "ConeInjection.H"
#include "ConeInjectionMP.H"
#include "FieldActivatedInjection.H"
#include "ManualInjection.H"
#include "NoInjection.H"
#include "PatchInjection.H"
#include "ReactingMultiphaseLookupTableInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingMultiphaseInjectionModels(ParcelType)                     \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        constGasThermoPhysics                                                 \
    );                                                                        \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        gasThermoPhysics                                                      \
    );                                                                        \
                                                                              \
    makeReactingMultiphaseInjectionModelThermoType                            \
    (                                                                         \
        ParcelType,                                                           \
        icoPoly8ThermoPhysics                                                 \
    );


#define makeReactingMultiphaseInjectionModelThermoType(ParcelType, ThermoType)\
                                                                              \
    makeInjectionModel(KinematicCloud<ParcelType<ThermoType> >);              \
                                                                              \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ConeInjection,                                                        \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ConeInjectionMP,                                                      \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        FieldActivatedInjection,                                              \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ManualInjection,                                                      \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        NoInjection,                                                          \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        PatchInjection,                                                       \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );                                                                        \
    makeInjectionModelThermoType                                              \
    (                                                                         \
        ReactingMultiphaseLookupTableInjection,                               \
        KinematicCloud,                                                       \
        ParcelType,                                                           \
        ThermoType                                                            \
    );


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
