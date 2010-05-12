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

#include "basicReactingParcel.H"
#include "ReactingCloud.H"

namespace Foam
{
    defineTemplateTypeNameAndDebug(Cloud<basicReactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicParcel<basicReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicParcel<basicReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ThermoParcel<basicReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ThermoParcel<basicReactingParcel>, 0);
    defineParcelTypeNameAndDebug(ReactingParcel<basicReactingParcel>, 0);
    defineTemplateTypeNameAndDebug(ReactingParcel<basicReactingParcel>, 0);

    defineParcelTypeNameAndDebug(KinematicCloud<basicReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(KinematicCloud<basicReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ThermoCloud<basicReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ThermoCloud<basicReactingParcel>, 0);

    defineParcelTypeNameAndDebug(ReactingCloud<basicReactingParcel>, 0);
//    defineTemplateTypeNameAndDebug(ReactingCloud<basicReactingParcel>, 0);

};


// ************************************************************************* //
