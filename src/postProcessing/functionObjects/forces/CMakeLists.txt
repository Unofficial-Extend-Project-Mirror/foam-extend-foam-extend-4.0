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
  forces/forces.C
  forces/forcesFunctionObject.C
  forceCoeffs/forceCoeffs.C
  forceCoeffs/forceCoeffsFunctionObject.C
)

set(sDoFRBM pointPatchFields/derived/sixDoFRigidBodyMotion)

list(APPEND SOURCES
  ${sDoFRBM}/sixDoFRigidBodyMotion.C
  ${sDoFRBM}/sixDoFRigidBodyMotionIO.C
  ${sDoFRBM}/sixDoFRigidBodyMotionState.C
  ${sDoFRBM}/sixDoFRigidBodyMotionStateIO.C
)

set(sDoFRBMR ${sDoFRBM}/sixDoFRigidBodyMotionRestraint)

list(APPEND SOURCES
  ${sDoFRBMR}/sixDoFRigidBodyMotionRestraint/sixDoFRigidBodyMotionRestraint.C
  ${sDoFRBMR}/sixDoFRigidBodyMotionRestraint/newSixDoFRigidBodyMotionRestraint.C
  ${sDoFRBMR}/linearAxialAngularSpring/linearAxialAngularSpring.C
  ${sDoFRBMR}/linearSpring/linearSpring.C
  ${sDoFRBMR}/sphericalAngularSpring/sphericalAngularSpring.C
  ${sDoFRBMR}/tabulatedAxialAngularSpring/tabulatedAxialAngularSpring.C
)

set(sDoFRBMC ${sDoFRBM}/sixDoFRigidBodyMotionConstraint)

list(APPEND SOURCES
  ${sDoFRBMC}/sixDoFRigidBodyMotionConstraint/sixDoFRigidBodyMotionConstraint.C
  ${sDoFRBMC}/sixDoFRigidBodyMotionConstraint/newSixDoFRigidBodyMotionConstraint.C
  ${sDoFRBMC}/fixedAxis/fixedAxis.C
  ${sDoFRBMC}/fixedLine/fixedLine.C
  ${sDoFRBMC}/fixedOrientation/fixedOrientation.C
  ${sDoFRBMC}/fixedPlane/fixedPlane.C
  ${sDoFRBMC}/fixedPoint/fixedPoint.C
  pointPatchFields/derived/sixDoFRigidBodyDisplacement/sixDoFRigidBodyDisplacementPointPatchVectorField.C
  pointPatchFields/derived/uncoupledSixDoFRigidBodyDisplacement/uncoupledSixDoFRigidBodyDisplacementPointPatchVectorField.C
)

add_foam_library(forces SHARED ${SOURCES})

target_link_libraries(forces PUBLIC incompressibleRASModels incompressibleLESModels compressibleRASModels compressibleLESModels)
