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
  engineTime/engineTime.C
  ignition/ignition.C
  ignition/ignitionIO.C
  ignition/ignitionSite.C
  ignition/ignitionSiteIO.C
  engineValve/engineValve.C
  engineVerticalValve/engineVerticalValve.C
  thoboisSlidingValve/thoboisSlidingValve.C
  dieselEngineValve/dieselEngineValve.C
  thoboisValve/thoboisValve.C
  accordionValve/accordionValve.C
  simpleEnginePiston/simpleEnginePiston.C
  enginePiston/enginePiston.C
  engineTopoChangerMesh/regionSide/regionSide.C
  engineTopoChangerMesh/attachDetachFunctions/attachDetachFunctions.C
  engineMesh/engineMesh/engineMesh.C
  engineMesh/engineMesh/newEngineMesh.C
  engineMesh/staticEngineMesh/staticEngineMesh.C
  engineMesh/layeredEngineMesh/layeredEngineMesh.C
  engineMesh/fvMotionSolverEngineMesh/fvMotionSolverEngineMesh.C
  engineMesh/tetDecompositionEngineMesh/tetDecompositionEngineMesh.C
  simpleEngineTopoFvMesh/simpleEngineTopoFvMesh.C
  simpleEngineTopoFvMesh/addSimpleEngineTopoFvMeshModifiers.C
  engineTopoChangerMesh/engineTopoChangerMesh/engineTopoChangerMesh.C
  engineTopoChangerMesh/engineTopoChangerMesh/newEngineTopoChangerMesh.C
  engineTopoChangerMesh/accordionEngineMesh/accordionEngineMesh.C
  engineTopoChangerMesh/accordionEngineMesh/accordionEngineMeshInitialize.C
  engineTopoChangerMesh/accordionEngineMesh/accordionEngineMeshMove.C
  engineTopoChangerMesh/accordionEngineMesh/addAccordionEngineMeshZones.C
  engineTopoChangerMesh/deformingEngineMesh/deformingEngineMesh.C
  engineTopoChangerMesh/deformingEngineMesh/deformingEngineMeshInitialize.C
  engineTopoChangerMesh/deformingEngineMesh/deformingEngineMeshMove.C
  engineTopoChangerMesh/deformingEngineMesh/addDeformingEngineMeshZones.C
  engineTopoChangerMesh/engineValveSliding/addEngineValveSlidingMeshModifiers.C
  engineTopoChangerMesh/engineValveSliding/engineValveSliding.C
  engineTopoChangerMesh/engineValveSliding/engineValveSlidingInitialize.C
  engineTopoChangerMesh/engineValveSliding/engineValveSlidingMove.C
  engineTopoChangerMesh/layerAR/layerAR.C
  engineTopoChangerMesh/layerAR/addLayerARMeshModifiers.C
  engineTopoChangerMesh/layerARGambit/layerARGambit.C
  engineTopoChangerMesh/layerARGambit/addLayerARGambitMeshModifiers.C
  engineTopoChangerMesh/pistonLayer/pistonLayer.C
  engineTopoChangerMesh/pistonLayer/addPistonLayerMeshModifiers.C
  engineTopoChangerMesh/pistonDeform/pistonDeform.C
  engineTopoChangerMesh/thoboisMesh/addThoboisMeshModifiers.C
  engineTopoChangerMesh/thoboisMesh/thoboisMesh.C
  engineTopoChangerMesh/thoboisMesh/thoboisMeshInitialize.C
  engineTopoChangerMesh/thoboisMesh/thoboisMeshMove.C
  engineTopoChangerMesh/twoStrokeEngine/twoStrokeEngine.C
  engineTopoChangerMesh/twoStrokeEngine/addTwoStrokeEngineModifiers.C
  engineTopoChangerMesh/twoStrokeEngine/twoStrokeEngineCalculate.C
  engineTopoChangerMesh/twoStrokeEngine/twoStrokeEngineInitialize.C
  engineTopoChangerMesh/twoStrokeEngine/twoStrokeEngineMove.C
  engineTopoChangerMesh/twoStrokeEngine/twoStrokeEngineCheckMotionFluxes.C
  engineTopoChangerMesh/simpleTwoStroke/simpleTwoStroke.C
  engineTopoChangerMesh/simpleTwoStroke/addSimpleTwoStrokeModifiers.C
  engineTopoChangerMesh/simpleTwoStroke/simpleTwoStrokeCalculate.C
  engineTopoChangerMesh/simpleTwoStroke/simpleTwoStrokeInitialize.C
  engineTopoChangerMesh/simpleTwoStroke/simpleTwoStrokeMove.C
  engineTopoChangerMesh/verticalValves/verticalValves.C
  engineTopoChangerMesh/verticalValves/addVerticalValvesMeshModifiers.C
  engineTopoChangerMesh/verticalValves/verticalValvesCalculate.C
  engineTopoChangerMesh/verticalValves/verticalValvesInitialize.C
  engineTopoChangerMesh/verticalValves/verticalValvesMove.C
  engineTopoChangerMesh/thoboisSliding/addThoboisSlidingMeshModifiers.C
  engineTopoChangerMesh/thoboisSliding/thoboisSliding.C
  engineTopoChangerMesh/thoboisSliding/thoboisSlidingCalculate.C
  engineTopoChangerMesh/thoboisSliding/thoboisSlidingInitialize.C
  engineTopoChangerMesh/thoboisSliding/thoboisSlidingMove.C
  engineTopoChangerMesh/pistonSliding/addPistonSlidingMeshModifiers.C
  engineTopoChangerMesh/pistonSliding/pistonSliding.C
  engineTopoChangerMesh/pistonSliding/pistonSlidingCalculate.C
  engineTopoChangerMesh/pistonSliding/pistonSlidingInitialize.C
  engineTopoChangerMesh/pistonSliding/pistonSlidingMove.C
  engineTopoChangerMesh/pistonRefine/pistonRefine.C
  engineTopoChangerMesh/noEngineMesh/noEngineMesh.C
)

#  derivedFvPatchFields/engineTimeVaryingUniformFixedValue/engineTimeVaryingUniformFixedValueFvPatchField.C
list(APPEND SOURCES
  derivedFvPatchFields/engineTimeVaryingUniformFixedValue/engineTimeVaryingUniformFixedValueFvPatchFields.C
  derivedFvPatchFields/engineMassFlowRateInletVelocity/engineMassFlowRateInletVelocityFvPatchVectorField.C
  derivedFvPatchFields/engineTimeVaryingTotalPressureFvPatchScalarField/engineTimeVaryingTotalPressureFvPatchScalarField.C
)

add_foam_library(engine SHARED ${SOURCES})

target_link_libraries(engine PUBLIC dynamicFvMesh fvMotionSolver)
