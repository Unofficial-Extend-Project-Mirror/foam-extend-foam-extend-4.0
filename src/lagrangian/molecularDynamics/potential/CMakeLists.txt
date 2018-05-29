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

set(potential potential)

list(APPEND SOURCES
  ${potential}/potential.C
)

set(pairPotential pairPotential)

list(APPEND SOURCES
  ${pairPotential}/pairPotentialList/pairPotentialList.C
  ${pairPotential}/basic/pairPotential.C
  ${pairPotential}/basic/pairPotentialIO.C
  ${pairPotential}/basic/newPairPotential.C
  ${pairPotential}/derived/lennardJones/lennardJones.C
  ${pairPotential}/derived/maitlandSmith/maitlandSmith.C
  ${pairPotential}/derived/azizChen/azizChen.C
  ${pairPotential}/derived/exponentialRepulsion/exponentialRepulsion.C
  ${pairPotential}/derived/coulomb/coulomb.C
  ${pairPotential}/derived/dampedCoulomb/dampedCoulomb.C
  ${pairPotential}/derived/noInteraction/noInteraction.C
)

set(energyScalingFunction energyScalingFunction)

list(APPEND SOURCES
  ${energyScalingFunction}/basic/energyScalingFunction.C
  ${energyScalingFunction}/basic/newEnergyScalingFunction.C
  ${energyScalingFunction}/derived/shifted/shifted.C
  ${energyScalingFunction}/derived/shiftedForce/shiftedForce.C
  ${energyScalingFunction}/derived/noScaling/noScaling.C
  ${energyScalingFunction}/derived/sigmoid/sigmoid.C
  ${energyScalingFunction}/derived/doubleSigmoid/doubleSigmoid.C
)

set(tetherPotential tetherPotential)

list(APPEND SOURCES
  ${tetherPotential}/tetherPotentialList/tetherPotentialList.C
  ${tetherPotential}/basic/tetherPotential.C
  ${tetherPotential}/basic/newTetherPotential.C
  ${tetherPotential}/derived/harmonicSpring/harmonicSpring.C
  ${tetherPotential}/derived/restrainedHarmonicSpring/restrainedHarmonicSpring.C
  ${tetherPotential}/derived/pitchForkRing/pitchForkRing.C
)

set(electrostaticPotential electrostaticPotential)

list(APPEND SOURCES
  ${electrostaticPotential}/electrostaticPotential.C
)

add_foam_library(potential SHARED ${SOURCES})

target_link_libraries(potential PUBLIC finiteVolume)
