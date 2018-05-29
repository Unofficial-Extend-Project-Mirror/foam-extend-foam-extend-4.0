# --------------------------------------------------------------------------
#   ========                 |
#   \      /  F ield         | foam-extend: Open Source CFD
#    \    /   O peration     | Version:     4.1
#     \  /    A nd           | Web:         http://www.foam-extend.org
#      \/     M anipulation  | For copyright notice see file Copyright
# --------------------------------------------------------------------------
# License
#     This file is part of foam-extend.
#
#     foam-extend is free software: you can redistribute it and/or modify it
#     under the terms of the GNU General Public License as published by the
#     Free Software Foundation, either version 3 of the License, or (at your
#     option) any later version.
#
#     foam-extend is distributed in the hope that it will be useful, but
#     WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#     General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.
#
# Description
#     CMakeLists.txt file for libraries and applications
#
# Author
#     Henrik Rusche, Wikki GmbH, 2017. All rights reserved
#
#
# --------------------------------------------------------------------------

set(parcel parcel)
set(spray spray)
set(injector injector)

set(atomizationModels spraySubModels/atomizationModel)
set(breakupModels spraySubModels/breakupModel)
set(dragModels spraySubModels/dragModel)
set(evaporationModels spraySubModels/evaporationModel)
set(heatTransferModels spraySubModels/heatTransferModel)
set(wallModels spraySubModels/wallModel)
set(collisionModels spraySubModels/collisionModel)
set(dispersionModels spraySubModels/dispersionModel)
set(injectorModels spraySubModels/injectorModel)

list(APPEND SOURCES
  ${parcel}/parcel.C
  ${parcel}/parcelFunctions.C
  ${parcel}/parcelIO.C
  ${parcel}/setRelaxationTimes.C
  ${spray}/spray.C
  ${spray}/sprayOps.C
  ${spray}/sprayInject.C
  ${spray}/sprayFunctions.C
  ${injector}/injector/injector.C
  ${injector}/injector/injectorIO.C
  ${injector}/injectorType/injectorType.C
  ${injector}/unitInjector/unitInjector.C
  ${injector}/multiHoleInjector/multiHoleInjector.C
  ${injector}/commonRailInjector/commonRailInjector.C
  ${injector}/swirlInjector/swirlInjector.C
  ${injector}/definedInjector/definedInjector.C
  ${atomizationModels}/atomizationModel/atomizationModel.C
  ${atomizationModels}/atomizationModel/newAtomizationModel.C
  ${atomizationModels}/LISA/LISA.C
  ${atomizationModels}/noAtomization/noAtomization.C
  ${atomizationModels}/blobsSheetAtomization/blobsSheetAtomization.C
  ${breakupModels}/breakupModel/newBreakupModel.C
  ${breakupModels}/breakupModel/breakupModel.C
  ${breakupModels}/noBreakup/noBreakup.C
  ${breakupModels}/reitzDiwakar/reitzDiwakar.C
  ${breakupModels}/reitzKHRT/reitzKHRT.C
  ${breakupModels}/SHF/SHF.C
  ${breakupModels}/TAB/TAB.C
  ${breakupModels}/ETAB/ETAB.C
  ${dragModels}/dragModel/newDragModel.C
  ${dragModels}/dragModel/dragModel.C
  ${dragModels}/noDragModel/noDragModel.C
  ${dragModels}/standardDragModel/standardDragModel.C
  ${evaporationModels}/evaporationModel/newEvaporationModel.C
  ${evaporationModels}/evaporationModel/evaporationModel.C
  ${evaporationModels}/noEvaporation/noEvaporation.C
  ${evaporationModels}/RutlandFlashBoil/RutlandFlashBoil.C
  ${evaporationModels}/standardEvaporationModel/standardEvaporationModel.C
  ${evaporationModels}/saturateEvaporationModel/saturateEvaporationModel.C
  ${heatTransferModels}/heatTransferModel/newHeatTransferModel.C
  ${heatTransferModels}/heatTransferModel/heatTransferModel.C
  ${heatTransferModels}/noHeatTransfer/noHeatTransfer.C
  ${heatTransferModels}/RanzMarshall/RanzMarshall.C
  ${injectorModels}/injectorModel/newInjectorModel.C
  ${injectorModels}/injectorModel/injectorModel.C
  ${injectorModels}/constant/constInjector.C
  ${injectorModels}/Chomiak/Chomiak.C
  ${injectorModels}/hollowCone/hollowCone.C
  ${injectorModels}/pressureSwirl/pressureSwirlInjector.C
  ${injectorModels}/definedHollowCone/definedHollowCone.C
  ${injectorModels}/definedPressureSwirl/definedPressureSwirl.C
  ${injectorModels}/blobsSwirl/blobsSwirlInjector.C
  ${wallModels}/wallModel/newWallModel.C
  ${wallModels}/wallModel/wallModel.C
  ${wallModels}/removeParcel/removeParcel.C
  ${wallModels}/reflectParcel/reflectParcel.C
  ${collisionModels}/collisionModel/collisionModel.C
  ${collisionModels}/collisionModel/newCollisionModel.C
  ${collisionModels}/noCollision/noCollision.C
  ${collisionModels}/ORourke/ORourkeCollisionModel.C
  ${collisionModels}/trajectoryModel/trajectoryModel.C
  ${dispersionModels}/dispersionModel/dispersionModel.C
  ${dispersionModels}/dispersionModel/newDispersionModel.C
  ${dispersionModels}/dispersionRASModel/dispersionRASModel.C
  ${dispersionModels}/dispersionLESModel/dispersionLESModel.C
  ${dispersionModels}/noDispersion/noDispersion.C
  ${dispersionModels}/gradientDispersionRAS/gradientDispersionRAS.C
  ${dispersionModels}/stochasticDispersionRAS/stochasticDispersionRAS.C
)

add_foam_library(dieselSpray SHARED ${SOURCES})

target_link_libraries(dieselSpray PUBLIC lagrangianBasic lagrangianIntermediate basicThermophysicalModels)
