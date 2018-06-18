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
  oversetRegion/oversetRegion.C

  oversetFringe/oversetFringe/oversetFringe.C
  oversetFringe/oversetFringe/newOversetFringe.C
  oversetFringe/manualFringe/manualFringe.C
  oversetFringe/faceCellsFringe/faceCellsFringe.C
  oversetFringe/overlapFringe/overlapFringe.C
  oversetFringe/compositeFringe/compositeFringe.C

  oversetFringe/overlapFringe/donorSuitability/donorSuitability/donorSuitability.C
  oversetFringe/overlapFringe/donorSuitability/donorSuitability/newDonorSuitability.C
  oversetFringe/overlapFringe/donorSuitability/noSuitability/noSuitability.C
  oversetFringe/overlapFringe/donorSuitability/patchDistance/patchDistance.C
  oversetFringe/overlapFringe/donorSuitability/cellVolumes/cellVolumes.C
  oversetFringe/overlapFringe/donorSuitability/faceArea/faceArea.C
  oversetFringe/overlapFringe/donorSuitability/cellBoundingBoxDiagonal/cellBoundingBoxDiagonal.C

  oversetInterpolation/oversetInterpolation/oversetInterpolation.C
  oversetInterpolation/oversetInterpolation/newOversetInterpolation.C
  oversetInterpolation/injectionInterpolation/injectionInterpolation.C
  oversetInterpolation/inverseDistanceInterpolation/inverseDistanceInterpolation.C
  oversetInterpolation/averageValueInterpolation/averageValueInterpolation.C

  oversetMesh/oversetMesh.C
  oversetMesh/oversetMeshAddressing.C

  oversetLduInterface/oversetLduInterface.C
  oversetLduInterfaceField/oversetLduInterfaceField.C

  oversetPolyPatch/oversetPolyPatch.C
  oversetPointPatch/oversetPointPatch.C

  oversetFvPatch/oversetFvPatch.C

  oversetFvPatchField/oversetFvPatchFields.C
  oversetFvPatchField/emptyOversetFvPatchField/emptyOversetFvPatchFields.C
  oversetFvsPatchField/oversetFvsPatchFields.C

  oversetAdjustPhi/oversetAdjustPhi.C
  oversetAdjustPhi/globalOversetAdjustPhi.C
  oversetAdjustPhi/regionWiseOversetAdjustPhi.C
)

add_foam_library(oversetMesh SHARED ${SOURCES})

target_link_libraries(oversetMesh PUBLIC sampling dynamicMesh)
