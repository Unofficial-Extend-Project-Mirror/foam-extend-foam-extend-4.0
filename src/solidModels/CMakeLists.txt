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

set(arbitraryCrack arbitraryCrack)
list(APPEND SOURCES
  ${arbitraryCrack}/faceCracker/faceCracker.C
  ${arbitraryCrack}/faceCracker/detachFaceCracker.C
  ${arbitraryCrack}/crackerFvMesh/crackerFvMesh.C
  ${arbitraryCrack}/solidCohesive/solidCohesiveFvPatchVectorField.C
  ${arbitraryCrack}/solidCohesiveFixedModeMix/solidCohesiveFixedModeMixFvPatchVectorField.C
  ${arbitraryCrack}/cohesive/cohesivePolyPatch.C
)

set(boundaryConditions fvPatchFields)
list(APPEND SOURCES
  ${boundaryConditions}/analyticalPlateHoleTraction/analyticalPlateHoleTractionFvPatchVectorField.C
  ${boundaryConditions}/fixedDisplacement/fixedDisplacementFvPatchVectorField.C
  ${boundaryConditions}/fixedDisplacementFixedRotation/fixedDisplacementFixedRotationFvPatchVectorField.C
  ${boundaryConditions}/fixedDisplacementZeroShear/fixedDisplacementZeroShearFvPatchVectorField.C
  ${boundaryConditions}/fixedRotation/fixedRotationFvPatchVectorField.C
  ${boundaryConditions}/solidContact/solidContactFvPatchVectorField.C
  ${boundaryConditions}/solidDirectionMixed/solidDirectionMixedFvPatchVectorField.C
  ${boundaryConditions}/solidTraction/solidTractionFvPatchVectorField.C
  ${boundaryConditions}/solidTractionFree/solidTractionFreeFvPatchVectorField.C
  ${boundaryConditions}/solidSymmetry/solidSymmetryFvPatchVectorField.C
  ${boundaryConditions}/timeVaryingSolidTraction/timeVaryingSolidTractionFvPatchVectorField.C
  ${boundaryConditions}/timeVaryingFixedDisplacement/timeVaryingFixedDisplacementFvPatchVectorField.C
  ${boundaryConditions}/timeVaryingFixedDisplacementZeroShear/timeVaryingFixedDisplacementZeroShearFvPatchVectorField.C
  ${boundaryConditions}/timeVaryingFixedRotation/timeVaryingFixedRotationFvPatchVectorField.C
  ${boundaryConditions}/transitionalParabolicVelocity/transitionalParabolicVelocityFvPatchVectorField.C
)

set(contactModels contactModels)
list(APPEND SOURCES
  ${contactModels}/normalContactModels/normalContactModel/normalContactModel.C
  ${contactModels}/normalContactModels/normalContactModel/newNormalContactModel.C
  ${contactModels}/normalContactModels/iterativePenalty/iterativePenalty.C
  ${contactModels}/normalContactModels/standardPenalty/standardPenalty.C
  ${contactModels}/normalContactModels/dirichletNeumann/dirichletNeumann.C
  ${contactModels}/frictionContactModels/frictionContactModel/frictionContactModel.C
  ${contactModels}/frictionContactModels/frictionContactModel/newFrictionContactModel.C
  ${contactModels}/frictionContactModels/frictionless/frictionless.C
  ${contactModels}/frictionContactModels/iterativePenaltyFriction/iterativePenaltyFriction.C
  ${contactModels}/frictionContactModels/standardPenaltyFriction/standardPenaltyFriction.C
  ${contactModels}/frictionContactModels/dirichletNeumannFriction/dirichletNeumannFriction.C
  ${contactModels}/frictionContactModels/frictionLaws/frictionLaw/frictionLaw.C
  ${contactModels}/frictionContactModels/frictionLaws/frictionLaw/newFrictionLaw.C
  ${contactModels}/frictionContactModels/frictionLaws/coulombFriction/coulombFriction.C
)

set(constitutiveModel constitutiveModel)
list(APPEND SOURCES
  ${constitutiveModel}/constitutiveModel.C
  ${constitutiveModel}/plasticityStressReturnMethods/plasticityStressReturn/plasticityStressReturn.C
  ${constitutiveModel}/plasticityStressReturnMethods/plasticityStressReturn/newPlasticityStressReturn.C
  ${constitutiveModel}/plasticityStressReturnMethods/aravasMises/aravasMises.C
  ${constitutiveModel}/plasticityStressReturnMethods/yamadaMises/yamadaMises.C
  ${constitutiveModel}/solidInterfaces/solidInterface/solidInterface.C
  ${constitutiveModel}/solidInterfaces/solidInterface/newSolidInterface.C
  ${constitutiveModel}/solidInterfaces/smallStrain/smallStrainSolidInterface.C
  ${constitutiveModel}/solidInterfaces/smallStrainCorrected/smallStrainCorrectedSolidInterface.C
  ${constitutiveModel}/solidInterfaces/smallStrainOrthotropic/smallStrainOrthotropicSolidInterface.C
  ${constitutiveModel}/solidInterfaces/smallStrainOrthotropicCorrected/smallStrainOrthotropicCorrectedSolidInterface.C
  ${constitutiveModel}/solidInterfaces/largeStrainUL/largeStrainULSolidInterface.C
  ${constitutiveModel}/solidInterfaces/largeStrainULCorrected/largeStrainULCorrectedSolidInterface.C
  ${constitutiveModel}/solidInterfaces/none/noneSolidInterface.C
  ${constitutiveModel}/tractionBoundaryGradient/tractionBoundaryGradient.C
)

set(cohesiveLaws ${constitutiveModel}/cohesiveLaws)
list(APPEND SOURCES
  ${cohesiveLaws}/cohesiveLaw/cohesiveLaw.C
  ${cohesiveLaws}/cohesiveLaw/newCohesiveLaw.C
  ${cohesiveLaws}/multiMaterialCohesiveLaw/multiMaterialCohesiveLaw.C
  ${cohesiveLaws}/dugdaleCohesiveLaw/dugdaleCohesiveLaw.C
  ${cohesiveLaws}/linearCohesiveLaw/linearCohesiveLaw.C
)

set(finiteVolume finiteVolume)
list(APPEND SOURCES
  ${finiteVolume}/gradSchemes/leastSquaresSolidInterfaceGrad/leastSquaresSolidInterfaceGrads.C
  ${finiteVolume}/gradSchemes/leastSquaresSolidInterfaceGrad/leastSquaresSolidInterfaceVectors.C
  nonLinearGeometry/nonLinearGeometry.C
)

set(rheologyLaws ${constitutiveModel}/rheologyLaws)
list(APPEND SOURCES
  ${rheologyLaws}/rheologyLaw/rheologyLaw.C
  ${rheologyLaws}/rheologyLaw/newRheologyLaw.C
  ${rheologyLaws}/linearElastic/linearElastic.C
  ${rheologyLaws}/elasticPlastic/elasticPlastic.C
  ${rheologyLaws}/linearElasticTabulatedPlastic/linearElasticTabulatedPlastic.C
  ${rheologyLaws}/multiMaterial/multiMaterial.C
  ${rheologyLaws}/orthotropicLinearElastic/orthotropicLinearElastic.C
  ${rheologyLaws}/PronyViscoelastic/PronyViscoelastic.C
)

add_foam_library(solidModels SHARED ${SOURCES})

target_link_libraries(solidModels PUBLIC topoChangerFvMesh finiteArea conjugateHeatTransfer)
