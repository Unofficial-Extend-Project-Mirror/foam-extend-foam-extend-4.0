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
  faMesh/faGlobalMeshData/faGlobalMeshData.C
  faMesh/faMesh.C
  faMesh/faMeshDemandDrivenData.C
  faMesh/faMeshUpdate.C
  faMesh/faBoundaryMesh/faBoundaryMesh.C
)

set(faPatches faMesh/faPatches)
list(APPEND SOURCES
  ${faPatches}/faPatch/faPatch.C
  ${faPatches}/faPatch/newFaPatch.C
  ${faPatches}/basic/coupled/coupledFaPatch.C
  ${faPatches}/constraint/empty/emptyFaPatch.C
  ${faPatches}/constraint/processor/processorFaPatch.C
  ${faPatches}/constraint/wedge/wedgeFaPatch.C
  ${faPatches}/constraint/cyclic/cyclicFaPatch.C
  ${faPatches}/constraint/symmetry/symmetryFaPatch.C
)

set(faMeshMapper faMesh/faMeshMapper)
list(APPEND SOURCES
  ${faMeshMapper}/faMeshMapper.C
  ${faMeshMapper}/faAreaMapper.C
  ${faMeshMapper}/faEdgeMapper.C
  ${faMeshMapper}/faPatchMapper.C
)

set(faPatchFields fields/faPatchFields)
list(APPEND SOURCES
  ${faPatchFields}/faPatchField/faPatchFields.C
  ${faPatchFields}/faPatchField/faPatchVectorNFields.C
)

set(basicFaPatchFields ${faPatchFields}/basic)
list(APPEND SOURCES
  ${basicFaPatchFields}/basicSymmetry/basicSymmetryFaPatchScalarField.C
  ${basicFaPatchFields}/calculated/calculatedFaPatchFields.C
  ${basicFaPatchFields}/calculated/calculatedFaPatchVectorNFields.C
  ${basicFaPatchFields}/coupled/coupledFaPatchFields.C
  ${basicFaPatchFields}/coupled/coupledFaPatchVectorNFields.C
  ${basicFaPatchFields}/zeroGradient/zeroGradientFaPatchFields.C
  ${basicFaPatchFields}/zeroGradient/zeroGradientFaPatchVectorNFields.C
  ${basicFaPatchFields}/fixedValue/fixedValueFaPatchFields.C
  ${basicFaPatchFields}/fixedValue/fixedValueFaPatchVectorNFields.C
  ${basicFaPatchFields}/fixedGradient/fixedGradientFaPatchFields.C
  ${basicFaPatchFields}/fixedGradient/fixedGradientFaPatchVectorNFields.C
  ${basicFaPatchFields}/mixed/mixedFaPatchFields.C
  ${basicFaPatchFields}/transform/transformFaPatchFields.C
  ${basicFaPatchFields}/transform/transformFaPatchScalarField.C
  ${basicFaPatchFields}/transform/transformFaPatchVectorNFields.C
)

set(constraintFaPatchFields ${faPatchFields}/constraint)
list(APPEND SOURCES
  ${constraintFaPatchFields}/empty/emptyFaPatchFields.C
  ${constraintFaPatchFields}/empty/emptyFaPatchVectorNFields.C
  ${constraintFaPatchFields}/processor/processorFaPatchFields.C
  ${constraintFaPatchFields}/processor/processorFaPatchScalarField.C
  ${constraintFaPatchFields}/processor/processorFaPatchVectorNFields.C
  ${constraintFaPatchFields}/wedge/wedgeFaPatchFields.C
  ${constraintFaPatchFields}/wedge/wedgeFaPatchScalarField.C
  ${constraintFaPatchFields}/wedge/wedgeFaPatchVectorNFields.C
  ${constraintFaPatchFields}/cyclic/cyclicFaPatchFields.C
  ${constraintFaPatchFields}/cyclic/cyclicFaPatchVectorNFields.C
  ${constraintFaPatchFields}/symmetry/symmetryFaPatchFields.C
)

set(derivedFaPatchFields ${faPatchFields}/derived)
list(APPEND SOURCES
  ${derivedFaPatchFields}/fixedValueOutflow/fixedValueOutflowFaPatchFields.C
  ${derivedFaPatchFields}/inletOutlet/inletOutletFaPatchFields.C
  ${derivedFaPatchFields}/slip/slipFaPatchFields.C
  ${derivedFaPatchFields}/edgeNormalFixedValue/edgeNormalFixedValueFaPatchVectorField.C
  ${derivedFaPatchFields}/timeVaryingUniformFixedValue/timeVaryingUniformFixedValueFaPatchFields.C
)

set(faePatchFields fields/faePatchFields)
list(APPEND SOURCES
  ${faePatchFields}/faePatchField/faePatchFields.C
)

set(basicFaePatchFields ${faePatchFields}/basic)
list(APPEND SOURCES
  ${basicFaePatchFields}/calculated/calculatedFaePatchFields.C
  ${basicFaePatchFields}/coupled/coupledFaePatchFields.C
  ${basicFaePatchFields}/fixedValue/fixedValueFaePatchFields.C
)

set(constraintFaePatchFields ${faePatchFields}/constraint)
list(APPEND SOURCES
  ${constraintFaePatchFields}/empty/emptyFaePatchFields.C
  ${constraintFaePatchFields}/processor/processorFaePatchFields.C
  ${constraintFaePatchFields}/wedge/wedgeFaePatchFields.C
  ${constraintFaePatchFields}/cyclic/cyclicFaePatchFields.C
  ${constraintFaePatchFields}/symmetry/symmetryFaePatchFields.C
)

set(faePatchVectorNFields ${faePatchFields}/faePatchVectorNFields)
list(APPEND SOURCES
  ${faePatchVectorNFields}/faePatchVectorNFields.C
  ${faePatchVectorNFields}/calculatedFaePatchVectorNFields.C
  ${faePatchVectorNFields}/emptyFaePatchVectorNFields.C
  ${faePatchVectorNFields}/wedgeFaePatchVectorNFields.C
  ${faePatchVectorNFields}/coupledFaePatchVectorNFields.C
  ${faePatchVectorNFields}/processorFaePatchVectorNFields.C

  fields/areaFields/areaFields.C
  fields/areaFields/areaVectorNFields.C
  fields/edgeFields/edgeFields.C
  fields/edgeFields/edgeVectorNFields.C
  faMatrices/faMatrices.C
  faMatrices/faScalarMatrix/faScalarMatrix.C
)

set(edgeInterpolation interpolation/edgeInterpolation)
list(APPEND SOURCES
  ${edgeInterpolation}/edgeInterpolation.C
  ${edgeInterpolation}/edgeInterpolationScheme/edgeInterpolationSchemes.C
)

set(schemes ${edgeInterpolation}/schemes)
list(APPEND SOURCES
  ${schemes}/linear/linearEdgeInterpolationMake.C
  ${schemes}/upwind/upwindEdgeInterpolationMake.C
  ${schemes}/linearUpwind/linearUpwindEdgeInterpolationMake.C
  ${schemes}/Gamma/GammaEdgeInterpolationMake.C
  ${schemes}/blended/blendedEdgeInterpolationMake.C
  finiteArea/fa/fa.C
  finiteArea/faSchemes/faSchemes.C
)

set(ddtSchemes finiteArea/ddtSchemes)
list(APPEND SOURCES
  ${ddtSchemes}/faDdtScheme/faDdtSchemes.C
  ${ddtSchemes}/steadyStateFaDdtScheme/steadyStateFaDdtSchemes.C
  ${ddtSchemes}/EulerFaDdtScheme/EulerFaDdtSchemes.C
  ${ddtSchemes}/backwardFaDdtScheme/backwardFaDdtSchemes.C
  ${ddtSchemes}/boundedBackwardFaDdtScheme/boundedBackwardFaDdtScheme.C
)

set(divSchemes finiteArea/divSchemes)
list(APPEND SOURCES
  finiteArea/fam/vectorFamDiv.C
  ${divSchemes}/faDivScheme/faDivSchemes.C
  ${divSchemes}/gaussFaDivScheme/gaussFaDivSchemes.C
)

set(gradSchemes finiteArea/gradSchemes)
list(APPEND SOURCES
  ${gradSchemes}/faGradScheme/faGradSchemes.C
  ${gradSchemes}/gaussFaGrad/gaussFaGrads.C
  ${gradSchemes}/leastSquaresFaGrad/leastSquaresFaVectors.C
  ${gradSchemes}/leastSquaresFaGrad/leastSquaresFaGrads.C
)

set(limitedGradSchemes ${gradSchemes}/limitedGradSchemes)
list(APPEND SOURCES
  ${limitedGradSchemes}/faceLimitedFaGrad/faceLimitedFaGrads.C
  ${limitedGradSchemes}/edgeLimitedFaGrad/edgeLimitedFaGrads.C
)

set(lnGradSchemes finiteArea/lnGradSchemes)
list(APPEND SOURCES
  ${lnGradSchemes}/lnGradScheme/lnGradSchemes.C
  ${lnGradSchemes}/correctedLnGrad/correctedLnGrads.C
  ${lnGradSchemes}/limitedLnGrad/limitedLnGrads.C
  ${lnGradSchemes}/fourthLnGrad/fourthLnGrads.C
)

set(laplacianSchemes finiteArea/laplacianSchemes)
list(APPEND SOURCES
  ${laplacianSchemes}/faLaplacianScheme/faLaplacianSchemes.C
  ${laplacianSchemes}/gaussFaLaplacianScheme/gaussFaLaplacianSchemes.C
)

set(convectionSchemes finiteArea/convectionSchemes)
list(APPEND SOURCES
  ${convectionSchemes}/faConvectionScheme/faConvectionSchemes.C
  ${convectionSchemes}/gaussFaConvectionScheme/gaussFaConvectionSchemes.C
)

set(general cfdTools/general)

set(solutionFaControl ${general}/solutionFaControl)
list(APPEND SOURCES
  ${solutionFaControl}/solutionFaControl/solutionFaControl.C
  ${solutionFaControl}/simpleFaControl/simpleFaControl.C
  ${solutionFaControl}/pimpleFaControl/pimpleFaControl.C
  ${solutionFaControl}/pisoFaControl/pisoFaControl.C
)

add_foam_library(finiteArea SHARED ${SOURCES})

target_link_libraries(finiteArea PUBLIC foam)
