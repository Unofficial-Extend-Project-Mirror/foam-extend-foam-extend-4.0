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
  fvMotionSolvers/fvMotionSolver/fvMotionSolver.C
  fvMotionSolvers/velocity/laplacian/velocityLaplacianFvMotionSolver.C
  fvMotionSolvers/velocity/refLaplacian/refVelocityLaplacianFvMotionSolver.C
  fvMotionSolvers/displacement/displacementFvMotionSolver/displacementFvMotionSolver.C
  fvMotionSolvers/displacement/interpolation/displacementInterpolationFvMotionSolver.C
  fvMotionSolvers/displacement/laplacian/displacementLaplacianFvMotionSolver.C
  fvMotionSolvers/displacement/SBRStress/displacementSBRStressFvMotionSolver.C
  fvMotionSolvers/velocity/componentLaplacian/velocityComponentLaplacianFvMotionSolver.C
  fvMotionSolvers/displacement/componentLaplacian/displacementComponentLaplacianFvMotionSolver.C
  motionDiffusivity/motionDiffusivity/motionDiffusivity.C
  motionDiffusivity/uniform/uniformDiffusivity.C
  motionDiffusivity/inverseDistance/inverseDistanceDiffusivity.C
  motionDiffusivity/inverseFaceDistance/inverseFaceDistanceDiffusivity.C
  motionDiffusivity/inversePointDistance/inversePointDistanceDiffusivity.C
  motionDiffusivity/inverseVolume/inverseVolumeDiffusivity.C
  motionDiffusivity/directional/directionalDiffusivity.C
  motionDiffusivity/motionDirectional/motionDirectionalDiffusivity.C
  motionDiffusivity/file/fileDiffusivity.C
  motionDiffusivity/manipulators/quadratic/quadraticDiffusivity.C
  motionDiffusivity/manipulators/exponential/exponentialDiffusivity.C
  fvPatchFields/derived/cellMotion/cellMotionFvPatchFields.C
  fvPatchFields/derived/surfaceSlipDisplacement/surfaceSlipDisplacementFvPatchFields.C
  pointPatchFields/derived/oscillatingVelocity/oscillatingVelocityPointPatchVectorField.C
  pointPatchFields/derived/angularOscillatingVelocity/angularOscillatingVelocityPointPatchVectorField.C
  pointPatchFields/derived/oscillatingDisplacement/oscillatingDisplacementPointPatchVectorField.C
  pointPatchFields/derived/angularOscillatingDisplacement/angularOscillatingDisplacementPointPatchVectorField.C
  pointPatchFields/derived/surfaceSlipDisplacement/surfaceSlipDisplacementPointPatchVectorField.C
  pointPatchFields/derived/surfaceDisplacement/surfaceDisplacementPointPatchVectorField.C
)

add_foam_library(fvMotionSolver SHARED ${SOURCES})

target_link_libraries(fvMotionSolver PUBLIC finiteVolume dynamicMesh)
