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
  fvMesh/fvMeshGeometry.C
  fvMesh/fvMesh.C
  fvMesh/singleCellFvMesh/singleCellFvMesh.C
  fvMesh/fvMeshSubset/fvMeshSubset.C
)

set(fvBoundaryMesh fvMesh/fvBoundaryMesh)
list(APPEND SOURCES
  ${fvBoundaryMesh}/fvBoundaryMesh.C
)

set(fvPatches fvMesh/fvPatches)
list(APPEND SOURCES
  ${fvPatches}/fvPatch/fvPatch.C
  ${fvPatches}/fvPatch/newFvPatch.C
)

set(basicFvPatches ${fvPatches}/basic)
list(APPEND SOURCES
  ${basicFvPatches}/coupled/coupledFvPatch.C
  ${basicFvPatches}/generic/genericFvPatch.C
)

set(constraintFvPatches ${fvPatches}/constraint)
list(APPEND SOURCES
  ${constraintFvPatches}/empty/emptyFvPatch.C
  ${constraintFvPatches}/symmetry/symmetryFvPatch.C
  ${constraintFvPatches}/wedge/wedgeFvPatch.C
  ${constraintFvPatches}/cyclic/cyclicFvPatch.C
  ${constraintFvPatches}/processor/processorFvPatch.C
  ${constraintFvPatches}/ggi/ggiFvPatch.C
  ${constraintFvPatches}/cyclicGgi/cyclicGgiFvPatch.C
  ${constraintFvPatches}/overlapGgi/overlapGgiFvPatch.C
  ${constraintFvPatches}/mixingPlane/mixingPlaneFvPatch.C
  ${constraintFvPatches}/regionCouple/regionCoupleFvPatch.C
)

set(derivedFvPatches ${fvPatches}/derived)
list(APPEND SOURCES
  ${derivedFvPatches}/wall/wallFvPatch.C
  ${derivedFvPatches}/directMapped/directMappedFvPatch.C
  ${derivedFvPatches}/cohesive/cohesiveFvPatch.C
  ${derivedFvPatches}/directMapped/directMappedWallFvPatch.C
)

set(wallDist fvMesh/wallDist)
list(APPEND SOURCES
  ${wallDist}/wallPointYPlus/wallPointYPlus.C
  ${wallDist}/nearWallDistNoSearch.C
  ${wallDist}/nearWallDist.C
  ${wallDist}/wallDist.C
  ${wallDist}/reflectionVectors.C
  ${wallDist}/wallDistReflection.C
)


set(fvMeshMapper fvMesh/fvMeshMapper)
list(APPEND SOURCES
  ${fvMeshMapper}/fvPatchMapper.C
  ${fvMeshMapper}/fvSurfaceMapper.C
)

set(extendedStencil fvMesh/extendedStencil)

set(cellToCell ${extendedStencil}/cellToCell)
list(APPEND SOURCES
  ${cellToCell}/fullStencils/cellToCellStencil.C
  ${cellToCell}/fullStencils/CFCCellToCellStencil.C
  ${cellToCell}/fullStencils/CPCCellToCellStencil.C
  ${cellToCell}/fullStencils/CECCellToCellStencil.C
)

set(cellToFace ${extendedStencil}/cellToFace)
list(APPEND SOURCES
  ${cellToFace}/fullStencils/cellToFaceStencil.C
  ${cellToFace}/fullStencils/CFCCellToFaceStencil.C
  ${cellToFace}/fullStencils/CECCellToFaceStencil.C
  ${cellToFace}/fullStencils/CPCCellToFaceStencil.C
  ${cellToFace}/fullStencils/FECCellToFaceStencil.C
  ${cellToFace}/extendedCellToFaceStencil.C
  ${cellToFace}/extendedCentredCellToFaceStencil.C
  ${cellToFace}/extendedUpwindCellToFaceStencil.C
  ${cellToFace}/MeshObjects/centredCECCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/centredCFCCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/centredCPCCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/centredFECCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/upwindCECCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/upwindCFCCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/upwindCPCCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/upwindFECCellToFaceStencilObject.C
  ${cellToFace}/MeshObjects/pureUpwindCFCCellToFaceStencilObject.C
)

set(faceToCell ${extendedStencil}/faceToCell)
list(APPEND SOURCES
  ${faceToCell}/fullStencils/faceToCellStencil.C
  ${faceToCell}/fullStencils/CFCFaceToCellStencil.C
  ${faceToCell}/extendedFaceToCellStencil.C
  ${faceToCell}/extendedCentredFaceToCellStencil.C
  ${faceToCell}/MeshObjects/centredCFCFaceToCellStencilObject.C
)


set(fvPatchFields fields/fvPatchFields)
list(APPEND SOURCES
  ${fvPatchFields}/fvPatchField/fvPatchFields.C
  ${fvPatchFields}/fvPatchField/fvPatchVectorNFields.C
)

set(basicFvPatchFields ${fvPatchFields}/basic)
list(APPEND SOURCES
#  ${basicFvPatchFields}/basicSymmetry/basicSymmetryFvPatchFields.C
  ${basicFvPatchFields}/basicSymmetry/basicSymmetryFvPatchScalarField.C
  ${basicFvPatchFields}/calculated/calculatedFvPatchFields.C
  ${basicFvPatchFields}/calculated/calculatedFvPatchVectorNFields.C
  ${basicFvPatchFields}/coupled/coupledFvPatchFields.C
  ${basicFvPatchFields}/coupled/coupledFvPatchVectorNFields.C
  ${basicFvPatchFields}/directionMixed/directionMixedFvPatchFields.C
  ${basicFvPatchFields}/fixedGradient/fixedGradientFvPatchFields.C
  ${basicFvPatchFields}/fixedGradient/fixedGradientFvPatchVectorNFields.C
  ${basicFvPatchFields}/fixedValue/fixedValueFvPatchFields.C
  ${basicFvPatchFields}/fixedValue/fixedValueFvPatchVectorNFields.C
  ${basicFvPatchFields}/generic/genericFvPatchFields.C
  ${basicFvPatchFields}/generic/genericFvPatchVectorNFields.C
  ${basicFvPatchFields}/mixed/mixedFvPatchFields.C
  ${basicFvPatchFields}/sliced/slicedFvPatchFields.C
  ${basicFvPatchFields}/transform/transformFvPatchFields.C
  ${basicFvPatchFields}/transform/transformFvPatchScalarField.C
  ${basicFvPatchFields}/transform/transformFvPatchVectorNFields.C
  ${basicFvPatchFields}/zeroGradient/zeroGradientFvPatchFields.C
  ${basicFvPatchFields}/zeroGradient/zeroGradientFvPatchVectorNFields.C
)

set(constraintFvPatchFields ${fvPatchFields}/constraint)
list(APPEND SOURCES
  ${constraintFvPatchFields}/cyclic/cyclicFvPatchFields.C
  ${constraintFvPatchFields}/cyclic/cyclicFvPatchVectorNFields.C
  ${constraintFvPatchFields}/empty/emptyFvPatchFields.C
  ${constraintFvPatchFields}/empty/emptyFvPatchVectorNFields.C
  ${constraintFvPatchFields}/jumpCyclic/jumpCyclicFvPatchFields.C
  ${constraintFvPatchFields}/processor/processorFvPatchFields.C
  ${constraintFvPatchFields}/processor/processorFvPatchScalarField.C
  ${constraintFvPatchFields}/processor/processorFvPatchVectorNFields.C
  ${constraintFvPatchFields}/symmetry/symmetryFvPatchFields.C
  ${constraintFvPatchFields}/wedge/wedgeFvPatchFields.C
  ${constraintFvPatchFields}/wedge/wedgeFvPatchScalarField.C
  ${constraintFvPatchFields}/wedge/wedgeFvPatchVectorNFields.C
  ${constraintFvPatchFields}/ggi/ggiFvPatchScalarField.C
  ${constraintFvPatchFields}/ggi/ggiFvPatchFields.C
  ${constraintFvPatchFields}/ggi/ggiFvPatchVectorNFields.C
  ${constraintFvPatchFields}/jumpGgi/jumpGgiFvPatchFields.C
  ${constraintFvPatchFields}/cyclicGgi/cyclicGgiFvPatchFields.C
  ${constraintFvPatchFields}/overlapGgi/overlapGgiFvPatchFields.C
  ${constraintFvPatchFields}/jumpOverlapGgi/jumpOverlapGgiFvPatchFields.C
  ${constraintFvPatchFields}/mixingPlane/mixingPlaneFvPatchFields.C
  ${constraintFvPatchFields}/jumpMixingPlane/jumpMixingPlaneFvPatchFields.C
  ${constraintFvPatchFields}/regionCoupling/regionCouplingFvPatchFields.C
)

set(derivedFvPatchFields ${fvPatchFields}/derived)
list(APPEND SOURCES
  ${derivedFvPatchFields}/activeBaffleVelocity/activeBaffleVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/advective/advectiveFvPatchFields.C
  ${derivedFvPatchFields}/directMappedFixedValue/directMappedFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/directMappedVelocityFluxFixedValue/directMappedVelocityFluxFixedValueFvPatchField.C
  ${derivedFvPatchFields}/fan/fanFvPatchFields.C
  ${derivedFvPatchFields}/buoyantPressure/buoyantPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/fixedFluxPressure/fixedFluxPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/fixedInternalValueFvPatchField/fixedInternalValueFvPatchFields.C
  ${derivedFvPatchFields}/fixedNormalSlip/fixedNormalSlipFvPatchFields.C
  ${derivedFvPatchFields}/fixedPressureCompressibleDensity/fixedPressureCompressibleDensityFvPatchScalarField.C
  ${derivedFvPatchFields}/fluxCorrectedVelocity/fluxCorrectedVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/flux/fluxFvPatchFields.C
  ${derivedFvPatchFields}/freestream/freestreamFvPatchFields.C
  ${derivedFvPatchFields}/freestreamPressure/freestreamPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/inletOutlet/inletOutletFvPatchFields.C
  ${derivedFvPatchFields}/inletOutletTotalTemperature/inletOutletTotalTemperatureFvPatchScalarField.C
  ${derivedFvPatchFields}/flowRateInletVelocity/flowRateInletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/movingWallVelocity/movingWallVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/translatingWallVelocity/translatingWallVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/rotatingWallVelocity/rotatingWallVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/movingRotatingWallVelocity/movingRotatingWallVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/oscillatingFixedValue/oscillatingFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/fixedMeanValue/fixedMeanValueFvPatchFields.C
  ${derivedFvPatchFields}/outletInlet/outletInletFvPatchFields.C
  ${derivedFvPatchFields}/partialSlip/partialSlipFvPatchFields.C
  ${derivedFvPatchFields}/pressureDirectedInletOutletVelocity/pressureDirectedInletOutletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/pressureNormalInletOutletVelocity/pressureNormalInletOutletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/pressureDirectedInletVelocity/pressureDirectedInletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/timeVaryingMappedPressureDirectedInletVelocity/timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/pressureInletOutletVelocity/pressureInletOutletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/pressureInletUniformVelocity/pressureInletUniformVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/pressureInletVelocity/pressureInletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/rotatingPressureInletOutletVelocity/rotatingPressureInletOutletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/rotatingTotalPressure/rotatingTotalPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/slip/slipFvPatchFields.C
  ${derivedFvPatchFields}/supersonicFreestream/supersonicFreestreamFvPatchVectorField.C
  ${derivedFvPatchFields}/surfaceNormalFixedValue/surfaceNormalFixedValueFvPatchVectorField.C
  ${derivedFvPatchFields}/syringePressure/syringePressureFvPatchScalarField.C
  ${derivedFvPatchFields}/parabolicVelocity/parabolicVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/timeVaryingMappedFixedValue/AverageIOFields.C
  ${derivedFvPatchFields}/timeVaryingMappedFixedValue/timeVaryingMappedFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/timeVaryingFlowRateInletVelocity/timeVaryingFlowRateInletVelocityFvPatchVectorField.C
  ${derivedFvPatchFields}/timeVaryingUniformFixedValue/timeVaryingUniformFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/timeVaryingUniformInletOutlet/timeVaryingUniformInletOutletFvPatchFields.C
  ${derivedFvPatchFields}/totalPressure/totalPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/timeVaryingMappedTotalPressure/timeVaryingMappedTotalPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/timeVaryingUniformTotalPressure/timeVaryingUniformTotalPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/totalTemperature/totalTemperatureFvPatchScalarField.C
  ${derivedFvPatchFields}/turbulentInlet/turbulentInletFvPatchFields.C
  ${derivedFvPatchFields}/turbulentIntensityKineticEnergyInlet/turbulentIntensityKineticEnergyInletFvPatchScalarField.C
  ${derivedFvPatchFields}/uniformFixedValue/uniformFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/waveTransmissive/waveTransmissiveFvPatchFields.C
  ${derivedFvPatchFields}/uniformDensityHydrostaticPressure/uniformDensityHydrostaticPressureFvPatchScalarField.C
  ${derivedFvPatchFields}/pulseFixedValue/pulseFixedValueFvPatchFields.C
  ${derivedFvPatchFields}/waveTransmissiveInlet/waveTransmissiveInletFvPatchFields.C
  ${derivedFvPatchFields}/noSlipWall/noSlipWallFvPatchVectorField.C
  ${derivedFvPatchFields}/noSlipMovingWall/noSlipMovingWallFvPatchVectorField.C
  ${derivedFvPatchFields}/blockSymmPlane/blockSymmPlaneFvPatchVectorField.C
)

set(fvsPatchFields fields/fvsPatchFields)
list(APPEND SOURCES
  ${fvsPatchFields}/fvsPatchField/fvsPatchFields.C
)

set(basicFvsPatchFields ${fvsPatchFields}/basic)
list(APPEND SOURCES
  ${basicFvsPatchFields}/calculated/calculatedFvsPatchFields.C
  ${basicFvsPatchFields}/coupled/coupledFvsPatchFields.C
  ${basicFvsPatchFields}/fixedValue/fixedValueFvsPatchFields.C
  ${basicFvsPatchFields}/sliced/slicedFvsPatchFields.C
)

set(constraintFvsPatchFields ${fvsPatchFields}/constraint)
list(APPEND SOURCES
  ${constraintFvsPatchFields}/cyclic/cyclicFvsPatchFields.C
  ${constraintFvsPatchFields}/empty/emptyFvsPatchFields.C
  ${constraintFvsPatchFields}/processor/processorFvsPatchFields.C
  ${constraintFvsPatchFields}/symmetry/symmetryFvsPatchFields.C
  ${constraintFvsPatchFields}/wedge/wedgeFvsPatchFields.C
  ${constraintFvsPatchFields}/ggi/ggiFvsPatchFields.C
  ${constraintFvsPatchFields}/cyclicGgi/cyclicGgiFvsPatchFields.C
  ${constraintFvsPatchFields}/overlapGgi/overlapGgiFvsPatchFields.C
  ${constraintFvsPatchFields}/mixingPlane/mixingPlaneFvsPatchFields.C
  ${constraintFvsPatchFields}/regionCoupling/regionCouplingFvsPatchFields.C
)

set(fvsPatchVectorNFields ${fvsPatchFields}/fvsPatchVectorNFields)
list(APPEND SOURCES
  ${fvsPatchVectorNFields}/fvsPatchVectorNFields.C
  ${fvsPatchVectorNFields}/calculatedFvsPatchVectorNFields.C
  ${fvsPatchVectorNFields}/emptyFvsPatchVectorNFields.C
  ${fvsPatchVectorNFields}/wedgeFvsPatchVectorNFields.C
  ${fvsPatchVectorNFields}/coupledFvsPatchVectorNFields.C
  ${fvsPatchVectorNFields}/processorFvsPatchVectorNFields.C
  fields/volFields/volFields.C
  fields/volFields/volVectorNFields.C
  fields/surfaceFields/surfaceFields.C
  fields/surfaceFields/surfaceVectorNFields.C
  fvMatrices/fvMatrices.C
  fvMatrices/fvScalarMatrix/fvScalarMatrix.C
  fvMatrices/solvers/MULES/MULES.C
  fvMatrices/solvers/GAMGSymSolver/GAMGAgglomerations/faceAreaPairGAMGAgglomeration/faceAreaPairGAMGAgglomeration.C
)

set(interpolation interpolation/interpolation)
list(APPEND SOURCES
  ${interpolation}/interpolation/interpolations.C
  ${interpolation}/interpolationCell/makeInterpolationCell.C
  ${interpolation}/interpolationCellPoint/cellPointWeight/cellPointWeight.C
  ${interpolation}/interpolationCellPoint/makeInterpolationCellPoint.C
  ${interpolation}/interpolationCellPointFace/makeInterpolationCellPointFace.C
  ${interpolation}/interpolationCellPointWallModified/cellPointWeightWallModified/cellPointWeightWallModified.C
  ${interpolation}/interpolationCellPointWallModified/makeInterpolationCellPointWallModified.C
  ${interpolation}/interpolationPointMVC/interpolationPointMVC.C
)

set(volPointInterpolation interpolation/volPointInterpolation)
list(APPEND SOURCES
  ${volPointInterpolation}/pointPatchInterpolation/pointPatchInterpolation.C
  ${volPointInterpolation}/volPointInterpolation.C
  ${volPointInterpolation}/leastSquaresVolPointInterpolation/leastSquaresVolPointInterpolation.C
)

set(pointVolInterpolation interpolation/pointVolInterpolation)
list(APPEND SOURCES
  ${pointVolInterpolation}/pointVolInterpolation.C
)

set(surfaceInterpolation interpolation/surfaceInterpolation)
list(APPEND SOURCES
  ${surfaceInterpolation}/surfaceInterpolation/surfaceInterpolation.C
  ${surfaceInterpolation}/surfaceInterpolationScheme/surfaceInterpolationSchemes.C
  ${surfaceInterpolation}/VectorNSurfaceInterpolationSchemes/VectorNSurfaceInterpolationSchemes.C
)

set(schemes ${surfaceInterpolation}/schemes)
list(APPEND SOURCES
  ${schemes}/linear/linear.C
  ${schemes}/midPoint/midPoint.C
  ${schemes}/downwind/downwind.C
  ${schemes}/weighted/weighted.C
  ${schemes}/cubic/cubic.C
  ${schemes}/skewCorrected/skewCorrectionVectors.C
  ${schemes}/skewCorrected/skewCorrected.C
  ${schemes}/leastSquaresSkewCorrected/leastSquaresSkewCorrected.C
  ${schemes}/outletStabilised/outletStabilised.C
  ${schemes}/reverseLinear/reverseLinear.C
  ${schemes}/clippedLinear/clippedLinear.C
  ${schemes}/harmonic/magLongDelta.C
  ${schemes}/harmonic/harmonic.C
  ${schemes}/fixedBlended/fixedBlended.C
  ${schemes}/localBlended/localBlended.C
  ${schemes}/localMax/localMax.C
  ${schemes}/localMin/localMin.C
  ${schemes}/linearFit/linearFit.C
  ${schemes}/biLinearFit/biLinearFit.C
  ${schemes}/quadraticLinearFit/quadraticLinearFit.C
  ${schemes}/quadraticFit/quadraticFit.C
  ${schemes}/quadraticLinearUpwindFit/quadraticLinearUpwindFit.C
  ${schemes}/quadraticUpwindFit/quadraticUpwindFit.C
  ${schemes}/cubicUpwindFit/cubicUpwindFit.C
)
#
# $(schemes)/quadraticLinearPureUpwindFit/quadraticLinearPureUpwindFit.C
list(APPEND SOURCES
  ${schemes}/linearPureUpwindFit/linearPureUpwindFit.C
)

set(limitedSchemes ${surfaceInterpolation}/limitedSchemes)
list(APPEND SOURCES
  ${limitedSchemes}/limitedSurfaceInterpolationScheme/limitedSurfaceInterpolationSchemes.C
  ${limitedSchemes}/upwind/upwind.C
  ${limitedSchemes}/blended/blended.C
  ${limitedSchemes}/linearUpwind/linearUpwind.C
  ${limitedSchemes}/linearUpwind/linearUpwindV.C
  ${limitedSchemes}/reconCentral/reconCentral.C
  ${limitedSchemes}/Gamma/Gamma.C
  ${limitedSchemes}/SFCD/SFCD.C
  ${limitedSchemes}/Minmod/Minmod.C
  ${limitedSchemes}/vanLeer/vanLeer.C
  ${limitedSchemes}/vanAlbada/vanAlbada.C
  ${limitedSchemes}/OSPRE/OSPRE.C
  ${limitedSchemes}/limitedLinear/limitedLinear.C
  ${limitedSchemes}/limitedCubic/limitedCubic.C
  ${limitedSchemes}/SuperBee/SuperBee.C
  ${limitedSchemes}/QUICK/QUICK.C
  ${limitedSchemes}/MUSCL/MUSCL.C
  ${limitedSchemes}/UMIST/UMIST.C
  ${limitedSchemes}/Phi/Phi.C
  ${limitedSchemes}/filteredLinear/filteredLinear.C
  ${limitedSchemes}/filteredLinear2/filteredLinear2.C
  ${limitedSchemes}/filteredLinear3/filteredLinear3.C
  ${limitedSchemes}/limitWith/limitWith.C
)

set(multivariateSchemes ${surfaceInterpolation}/multivariateSchemes)
list(APPEND SOURCES
  ${multivariateSchemes}/multivariateSurfaceInterpolationScheme/multivariateSurfaceInterpolationSchemes.C
  ${multivariateSchemes}/multivariateSelectionScheme/multivariateSelectionSchemes.C
  ${multivariateSchemes}/upwind/multivariateUpwind.C
  ${multivariateSchemes}/Gamma/multivariateGamma.C
  ${multivariateSchemes}/vanLeer/multivariateVanLeer.C
  ${multivariateSchemes}/Minmod/multivariateMinmod.C
  ${multivariateSchemes}/SuperBee/multivariateSuperBee.C
  ${multivariateSchemes}/MUSCL/multivariateMUSCL.C
  ${multivariateSchemes}/limitedLinear/multivariateLimitedLinear.C
  ${multivariateSchemes}/limitedCubic/multivariateLimitedCubic.C
  finiteVolume/fv/fv.C
  finiteVolume/fvSchemes/fvSchemes.C
)

set(ddtSchemes finiteVolume/ddtSchemes)
list(APPEND SOURCES
  ${ddtSchemes}/ddtScheme/ddtSchemes.C
  ${ddtSchemes}/steadyStateDdtScheme/steadyStateDdtSchemes.C
  ${ddtSchemes}/EulerDdtScheme/EulerDdtSchemes.C
  ${ddtSchemes}/CoEulerDdtScheme/CoEulerDdtSchemes.C
  ${ddtSchemes}/SLTSDdtScheme/SLTSDdtSchemes.C
  ${ddtSchemes}/backwardDdtScheme/backwardDdtSchemes.C
  ${ddtSchemes}/boundedBackwardDdtScheme/boundedBackwardDdtScheme.C
  ${ddtSchemes}/boundedBackwardDdtScheme/boundedBackwardDdtSchemes.C
  ${ddtSchemes}/CrankNicolsonDdtScheme/CrankNicolsonDdtSchemes.C
  ${ddtSchemes}/steadyInertialDdtScheme/steadyInertialDdtSchemes.C
)

set(d2dt2Schemes finiteVolume/d2dt2Schemes)
list(APPEND SOURCES
  ${d2dt2Schemes}/d2dt2Scheme/d2dt2Schemes.C
  ${d2dt2Schemes}/steadyStateD2dt2Scheme/steadyStateD2dt2Schemes.C
  ${d2dt2Schemes}/EulerD2dt2Scheme/EulerD2dt2Schemes.C
  ${d2dt2Schemes}/backwardD2dt2Scheme/backwardD2dt2Schemes.C
)

set(divSchemes finiteVolume/divSchemes)
list(APPEND SOURCES
  ${divSchemes}/divScheme/divSchemes.C
  ${divSchemes}/gaussDivScheme/gaussDivSchemes.C
  ${divSchemes}/gaussDivScheme/vectorGaussDivScheme.C
)

set(gradSchemes finiteVolume/gradSchemes)
list(APPEND SOURCES
  ${gradSchemes}/gradScheme/gradSchemes.C
  ${gradSchemes}/gaussGrad/scalarGaussGrad.C
  ${gradSchemes}/gaussGrad/gaussGrads.C
  ${gradSchemes}/beGaussGrad/beGaussGrads.C
  ${gradSchemes}/leastSquaresGrad/leastSquaresVectors.C
  ${gradSchemes}/leastSquaresGrad/scalarLeastSquaresGrad.C
  ${gradSchemes}/leastSquaresGrad/leastSquaresGrads.C
  ${gradSchemes}/extendedLeastSquaresGrad/extendedLeastSquaresVectors.C
  ${gradSchemes}/extendedLeastSquaresGrad/extendedLeastSquaresGrads.C
  ${gradSchemes}/fourthGrad/fourthGrads.C
)

set(limitedGradSchemes ${gradSchemes}/limitedGradSchemes)
list(APPEND SOURCES
  ${limitedGradSchemes}/faceLimitedGrad/faceLimitedGrads.C
  ${limitedGradSchemes}/cellLimitedGrad/cellLimitedGrads.C
  ${limitedGradSchemes}/faceMDLimitedGrad/faceMDLimitedGrads.C
  ${limitedGradSchemes}/cellMDLimitedGrad/cellMDLimitedGrads.C
)

set(snGradSchemes finiteVolume/snGradSchemes)
list(APPEND SOURCES
  ${snGradSchemes}/snGradScheme/snGradSchemes.C
  ${snGradSchemes}/correctedSnGrad/correctedSnGrads.C
  ${snGradSchemes}/limitedSnGrad/limitedSnGrads.C
  ${snGradSchemes}/uncorrectedSnGrad/uncorrectedSnGrads.C
  ${snGradSchemes}/orthogonalSnGrad/orthogonalSnGrads.C
  ${snGradSchemes}/skewCorrectedSnGrad/skewCorrectedSnGrads.C
  ${snGradSchemes}/fourthSnGrad/fourthSnGrads.C
)
#
# $(snGradSchemes)/quadraticFitSnGrad/quadraticFitSnGradData.C
# $(snGradSchemes)/quadraticFitSnGrad/quadraticFitSnGrads.C

set(convectionSchemes finiteVolume/convectionSchemes)
list(APPEND SOURCES
  ${convectionSchemes}/convectionScheme/convectionSchemes.C
  ${convectionSchemes}/gaussConvectionScheme/gaussConvectionSchemes.C
  ${convectionSchemes}/noConvectionScheme/noConvectionSchemes.C
  ${convectionSchemes}/explicitConvectionScheme/explicitConvectionSchemes.C
  ${convectionSchemes}/multivariateGaussConvectionScheme/multivariateGaussConvectionSchemes.C
)

set(adjConvectionSchemes finiteVolume/adjConvectionSchemes)
list(APPEND SOURCES
  ${adjConvectionSchemes}/adjConvectionScheme/adjConvectionSchemes.C
  ${adjConvectionSchemes}/noAdjConvectionScheme/noAdjConvectionSchemes.C
  ${adjConvectionSchemes}/explicitAdjConvectionScheme/explicitAdjConvectionSchemes.C
)

set(laplacianSchemes finiteVolume/laplacianSchemes)
list(APPEND SOURCES
  ${laplacianSchemes}/laplacianScheme/laplacianSchemes.C
  ${laplacianSchemes}/gaussLaplacianScheme/gaussLaplacianSchemes.C
  ${laplacianSchemes}/noLaplacianScheme/noLaplacianSchemes.C
  finiteVolume/fvc/fvcMeshPhi.C
)

set(general cfdTools/general)
list(APPEND SOURCES
  ${general}/findRefCell/findRefCell.C
  ${general}/adjustPhi/adjustPhi.C
  ${general}/bound/bound.C
)

set(solutionControl ${general}/solutionControl)
list(APPEND SOURCES
  ${solutionControl}/solutionControl/solutionControl.C
  ${solutionControl}/simpleControl/simpleControl.C
  ${solutionControl}/pimpleControl/pimpleControl.C
  ${solutionControl}/pisoControl/pisoControl.C
)

set(porousMedia ${general}/porousMedia)
list(APPEND SOURCES
  ${porousMedia}/porousZone.C
  ${porousMedia}/porousZones.C
)

set(MRF ${general}/MRF)
list(APPEND SOURCES
  ${MRF}/MRFZone.C
  ${MRF}/MRFZones.C
)

set(SRF ${general}/SRF)
list(APPEND SOURCES
  ${SRF}/SRFModel/SRFModel/SRFModel.C
  ${SRF}/SRFModel/SRFModel/newSRFModel.C
  ${SRF}/SRFModel/rpm/rpm.C
  ${SRF}/derivedFvPatchFields/SRFVelocityFvPatchVectorField/SRFVelocityFvPatchVectorField.C
  ${SRF}/derivedFvPatchFields/SRFSurfaceNormalVelocityFvPatchVectorField/SRFSurfaceNormalVelocityFvPatchVectorField.C
  ${SRF}/derivedFvPatchFields/SRFTotalTemperature/SRFTotalTemperatureFvPatchScalarField.C
  ${SRF}/derivedFvPatchFields/SRFFlowRateInletVelocity/SRFFlowRateInletVelocityFvPatchVectorField.C
  ${SRF}/derivedFvPatchFields/SRFTotalPressure/SRFTotalPressureFvPatchScalarField.C
)

set(fieldSources ${general}/fieldSources)
list(APPEND SOURCES
  ${fieldSources}/pressureGradientExplicitSource/pressureGradientExplicitSource.C
  ${fieldSources}/timeActivatedExplicitSource/timeActivatedExplicitSource.C
)

add_foam_library(finiteVolume SHARED ${SOURCES})

target_link_libraries(finiteVolume PUBLIC meshTools)

target_include_directories(finiteVolume PUBLIC
$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cfdTools>
)
