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

list(APPEND SOURCES
  codedFunctionObject/codedFunctionObject.C
  CourantNo/CourantNo.C
  CourantNo/CourantNoFunctionObject.C
  Lambda2/Lambda2.C
  Lambda2/Lambda2FunctionObject.C
  Peclet/Peclet.C
  Peclet/PecletFunctionObject.C
  Q/Q.C
  Q/QFunctionObject.C
  staticPressure/staticPressure.C
  staticPressure/staticPressureFunctionObject.C
  dsmcFields/dsmcFields.C
  dsmcFields/dsmcFieldsFunctionObject.C
  MachNumber/MachNumber.C
  divFlux/divFlux.C
  pressureTools/pressureTools.C
  pressureTools/pressureToolsFunctionObject.C
  residuals/residuals.C
  residuals/residualsFunctionObject.C
  scalarTransport/scalarTransport.C
  scalarTransport/scalarTransportFunctionObject.C
  timeActivatedFileUpdate/timeActivatedFileUpdate.C
  timeActivatedFileUpdate/timeActivatedFileUpdateFunctionObject.C
  turbulenceFields/turbulenceFields.C
  turbulenceFields/turbulenceFieldsFunctionObject.C
  vorticity/vorticity.C
  vorticity/vorticityFunctionObject.C
  wallShearStress/wallShearStress.C
  wallShearStress/wallShearStressFunctionObject.C
  setTimeStep/setTimeStepFunctionObject.C
)

add_foam_library(utilityFunctionObjects SHARED ${SOURCES})

target_link_libraries(utilityFunctionObjects PUBLIC basicThermophysicalModels sampling dsmc incompressibleTurbulenceModel compressibleTurbulenceModel)
