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

set(tetPolyMesh tetPolyMesh)
set(tetPolyBoundaryMesh ${tetPolyMesh}/tetPolyBoundaryMesh)
set(tetPolyPatches ${tetPolyMesh}/tetPolyPatches)
set(tetPolyPatch ${tetPolyPatches}/tetPolyPatch)
set(faceTetPolyPatch ${tetPolyPatches}/faceTetPolyPatch)
set(basicTetPolyPatches ${tetPolyPatches}/basic)
set(constraintTetPolyPatches ${tetPolyPatches}/constraint)
set(derivedTetPolyPatches ${tetPolyPatches}/derived)
set(coupledTetPolyPatch ${basicTetPolyPatches}/coupled)
set(coupledFaceTetPolyPatch ${basicTetPolyPatches}/coupledFace)
set(emptyTetPolyPatch ${constraintTetPolyPatches}/empty)
set(symmetryTetPolyPatch ${constraintTetPolyPatches}/symmetry)
set(wedgeTetPolyPatch ${constraintTetPolyPatches}/wedge)
set(processorTetPolyPatch ${constraintTetPolyPatches}/processor)
set(ggiTetPolyPatch ${constraintTetPolyPatches}/ggi)
set(cyclicGgiTetPolyPatch ${constraintTetPolyPatches}/cyclicGgi)
set(mixingPlaneTetPolyPatch ${constraintTetPolyPatches}/mixingPlane)
set(globalTetPolyPatch ${constraintTetPolyPatches}/global)
set(wallTetPolyPatch ${derivedTetPolyPatches}/wall)
set(directMappedTetPolyPatch ${derivedTetPolyPatches}/directMapped)

set(MapTetFemFields ${tetPolyMesh}/MapTetFemFields)

set(tetPolyPatchInterpolation tetPolyPatchInterpolation)

list(APPEND SOURCES
  ${tetPolyPatch}/tetPolyPatch.C
  ${faceTetPolyPatch}/faceTetPolyPatch.C
  ${faceTetPolyPatch}/newFaceTetPolyPatch.C
  ${coupledTetPolyPatch}/coupledTetPolyPatch.C
  ${emptyTetPolyPatch}/emptyTetPolyPatch.C
  ${symmetryTetPolyPatch}/symmetryTetPolyPatch.C
  ${wedgeTetPolyPatch}/wedgeTetPolyPatch.C
  ${coupledFaceTetPolyPatch}/coupledFaceTetPolyPatch.C
  ${processorTetPolyPatch}/processorTetPolyPatchM.C
  ${processorTetPolyPatch}/calcProcessorTetPolyPatchAddr.C
  ${processorTetPolyPatch}/calcProcessorTetPolyPatchPointAddr.C
  ${ggiTetPolyPatch}/ggiTetPolyPatch.C
  ${cyclicGgiTetPolyPatch}/cyclicGgiTetPolyPatch.C
  ${mixingPlaneTetPolyPatch}/mixingPlaneTetPolyPatch.C
  ${globalTetPolyPatch}/globalTetPolyPatch.C
  ${globalTetPolyPatch}/calcGlobalTetPolyPatchAddr.C
  ${wallTetPolyPatch}/wallTetPolyPatch.C
  ${directMappedTetPolyPatch}/directMappedTetPolyPatch.C
  ${directMappedTetPolyPatch}/directMappedWallTetPolyPatch.C
  ${tetPolyBoundaryMesh}/tetPolyBoundaryMesh.C
  ${tetPolyMesh}/tetPolyMeshLduAddressing.C
  ${tetPolyMesh}/tetPolyMesh.C
  ${tetPolyMesh}/calcTetPolyMeshGeometry.C
  ${tetPolyMesh}/calcTetPolyMeshAddressing.C
  ${tetPolyMesh}/addParallelPointPatch.C
  ${tetPolyMesh}/calcTetPolyMeshParPointData.C
  ${MapTetFemFields}/tetPointMapper.C
  ${MapTetFemFields}/tetPolyPatchMapper.C
  ${tetPolyPatchInterpolation}/tetPolyPatchInterpolation.C
)

set(tetPolyPatchFields fields/tetPolyPatchFields)
list(APPEND SOURCES
  ${tetPolyPatchFields}/tetPolyPatchField/tetPolyPatchFields.C
)

set(basicTetPolyPatchFields ${tetPolyPatchFields}/basic)
list(APPEND SOURCES
  ${basicTetPolyPatchFields}/calculated/calculatedTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/coupled/coupledTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/coupledFace/coupledFaceTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/generic/genericTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/value/valueTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/fixedValue/fixedValueTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/mixed/mixedTetPolyPatchFields.C
  ${basicTetPolyPatchFields}/zeroGradient/zeroGradientTetPolyPatchFields.C
)

set(constraintTetPolyPatchFields ${tetPolyPatchFields}/constraint)
list(APPEND SOURCES
  ${constraintTetPolyPatchFields}/empty/emptyTetPolyPatchFields.C
  ${constraintTetPolyPatchFields}/wedge/wedgeTetPolyPatchFields.C
  ${constraintTetPolyPatchFields}/symmetry/symmetryTetPolyPatchFields.C
  ${constraintTetPolyPatchFields}/processor/processorTetPolyPatchFields.C
)

set(derivedTetPolyPatchFields ${tetPolyPatchFields}/derived)
list(APPEND SOURCES
  ${derivedTetPolyPatchFields}/global/globalTetPolyPatchFields.C
  ${derivedTetPolyPatchFields}/componentMixed/componentMixedTetPolyPatchVectorField.C
  ${derivedTetPolyPatchFields}/oscillatingFixedValue/oscillatingFixedValueTetPolyPatchFields.C
  ${derivedTetPolyPatchFields}/timeVaryingUniformFixedValue/timeVaryingUniformFixedValueTetPolyPatchFields.C
  ${derivedTetPolyPatchFields}/slip/slipTetPolyPatchFields.C
)

set(elementPatchFields fields/elementPatchFields)
list(APPEND SOURCES
  ${elementPatchFields}/elementPatchField/elementPatchFields.C
)

set(basicElementPatchFields ${elementPatchFields}/basic)
list(APPEND SOURCES
  ${basicElementPatchFields}/calculated/calculatedElementPatchFields.C
  ${basicElementPatchFields}/coupled/coupledElementPatchFields.C
  ${basicElementPatchFields}/coupledFace/coupledFaceElementPatchFields.C
  ${basicElementPatchFields}/value/valueElementPatchFields.C
)

set(constraintElementPatchFields ${elementPatchFields}/constraint)
list(APPEND SOURCES
  ${constraintElementPatchFields}/empty/emptyElementPatchFields.C
  ${constraintElementPatchFields}/wedge/wedgeElementPatchFields.C
  ${constraintElementPatchFields}/symmetry/symmetryElementPatchFields.C
  ${constraintElementPatchFields}/processor/processorElementPatchFields.C
  ${constraintElementPatchFields}/global/globalElementPatchFields.C
  fields/tetPointFields/tetPointFields.C
  fields/elementFields/elementFields.C
  tetFemMatrix/tetFemScalarMatrix.C
  tetFemMatrix/tetFemMatrices.C
)

add_foam_library(tetFiniteElement SHARED ${SOURCES})

target_link_libraries(tetFiniteElement PUBLIC meshTools)
