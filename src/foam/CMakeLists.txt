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
#     Free Software Foundation, either version 3 of the License, or {at your
#     option} any later version.
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

list(APPEND sources
  global/global.C
  global/dimensionedConstants/dimensionedConstants.C
  global/argList/argList.C
  global/clock/clock.C
  global/controlSwitches/debugSwitch.C
  global/controlSwitches/infoSwitch.C
  global/controlSwitches/optimisationSwitch.C
  global/controlSwitches/tolerancesSwitch.C
  global/controlSwitches/constantsSwitch.C
  global/profiling/profilingInfo.C
  global/profiling/profilingPool.C
  global/profiling/profilingStack.C
  global/profiling/profilingTrigger.C
)

set(bools primitives/bools)
list(APPEND sources
  ${bools}/bool/bool.C
  ${bools}/bool/boolIO.C
  ${bools}/Switch/Switch.C
  ${bools}/Switch/SwitchIO.C
  primitives/char/charIO.C
)

set(ints primitives/ints)
list(APPEND sources
  ${ints}/int/intIO.C
  ${ints}/uint/uintIO.C
  ${ints}/long/longIO.C
  ${ints}/longLong/longLongIO.C
  ${ints}/ulong/ulongIO.C
  ${ints}/label/label.C
  ${ints}/uLabel/uLabel.C
  primitives/Scalar/doubleScalar/doubleScalar.C
  primitives/Scalar/longDoubleScalar/longDoubleScalar.C
  primitives/Scalar/floatScalar/floatScalar.C
  primitives/Scalar/scalar/scalar.C
  primitives/DiagTensor/diagTensor/diagTensor.C
  primitives/DiagTensor/labelDiagTensor/labelDiagTensor.C
  primitives/SphericalTensor/sphericalTensor/sphericalTensor.C
  primitives/SphericalTensor/labelSphericalTensor/labelSphericalTensor.C
  primitives/SymmTensor/labelSymmTensor/labelSymmTensor.C
  primitives/SymmTensor/symmTensor/symmTensor.C
  primitives/SymmTensor4thOrder/labelSymmTensor4thOrder/labelSymmTensor4thOrder.C
  primitives/SymmTensor4thOrder/symmTensor4thOrder/symmTensor4thOrder.C
  primitives/Tensor/labelTensor/labelTensor.C
  primitives/Tensor/tensor/tensor.C
  primitives/Vector/complexVector/complexVector.C
  primitives/Vector/labelVector/labelVector.C
  primitives/Vector/vector/vector.C
  primitives/Tensor2D/tensor2D/tensor2D.C
  primitives/SphericalTensor2D/sphericalTensor2D/sphericalTensor2D.C
  primitives/SymmTensor2D/symmTensor2D/symmTensor2D.C
  primitives/Vector2D/vector2D/vector2D.C
  primitives/complex/complex.C
  primitives/quaternion/quaternion.C
  primitives/septernion/septernion.C
)

set(strings primitives/strings)
list(APPEND sources
  ${strings}/string/foamString.C
  ${strings}/string/foamStringIO.C
  ${strings}/word/word.C
  ${strings}/word/wordIO.C
  ${strings}/fileName/fileName.C
  ${strings}/fileName/fileNameIO.C
  ${strings}/keyType/keyTypeIO.C
  ${strings}/wordRe/wordReIO.C
  primitives/hashes/Hasher/Hasher.C
)

set(sha1 primitives/hashes/SHA1)
list(APPEND sources
  ${sha1}/SHA1.C
  ${sha1}/SHA1Digest.C
  coordinateSystems/coordinateSystem.C
  coordinateSystems/coordinateSystems.C
  coordinateSystems/parabolicCylindricalCS.C
  coordinateSystems/ellipticCylindricalCS.C
  coordinateSystems/toroidalCS.C
  coordinateSystems/newCoordinateSystem.C
  coordinateSystems/cylindricalCS.C
  coordinateSystems/sphericalCS.C
  coordinateSystems/coordinateRotation/coordinateRotation.C
  coordinateSystems/coordinateRotation/EulerCoordinateRotation.C
  coordinateSystems/coordinateRotation/STARCDCoordinateRotation.C
  coordinateSystems/coordinateRotation/axisCoordinateRotation.C
  primitives/random/Random.C
  containers/HashTables/HashTable/HashTableCore.C
  containers/HashTables/StaticHashTable/StaticHashTableCore.C
  containers/Lists/SortableList/ParSortableListName.C
  containers/Lists/PackedList/PackedListName.C
  containers/Lists/ListOps/ListOps.C
  containers/LinkedLists/linkTypes/SLListBase/SLListBase.C
  containers/LinkedLists/linkTypes/DLListBase/DLListBase.C
)

set(primitiveLists primitives/Lists)
list(APPEND sources
  ${primitiveLists}/boolList.C
  ${primitiveLists}/diagTensorList.C
  ${primitiveLists}/labelIOList.C
  ${primitiveLists}/scalarList.C
  ${primitiveLists}/scalarIOList.C
  ${primitiveLists}/vectorList.C
  ${primitiveLists}/sphericalTensorList.C
  ${primitiveLists}/symmTensorList.C
  ${primitiveLists}/symmTensor4thOrderList.C
  ${primitiveLists}/tensorList.C
  ${primitiveLists}/VectorNLists.C
  ${primitiveLists}/wordIOList.C
)

set(Streams db/IOstreams)
list(APPEND sources
  ${Streams}/token/tokenIO.C
)

set(IOstreams ${Streams}/IOstreams)
list(APPEND sources
  ${IOstreams}/IOstream.C
  ${IOstreams}/Istream.C
  ${IOstreams}/Ostream.C
)

set(Sstreams ${Streams}/Sstreams)
list(APPEND sources
  ${Sstreams}/ISstream.C
  ${Sstreams}/OSstream.C
  ${Sstreams}/SstreamsPrint.C
  ${Sstreams}/readHexLabel.C
  ${Sstreams}/prefixOSstream.C
)

set(gzstream ${Streams}/gzstream)
list(APPEND sources
  ${gzstream}/gzstream.C
)

set(Fstreams ${Streams}/Fstreams)
list(APPEND sources
  ${Fstreams}/IFstream.C
  ${Fstreams}/OFstream.C
)

set(Tstreams ${Streams}/Tstreams)
list(APPEND sources
  ${Tstreams}/ITstream.C
)

set(StringStreams ${Streams}/StringStreams)
list(APPEND sources
  ${StringStreams}/StringStreamsPrint.C
)

set(Pstreams ${Streams}/Pstreams)
list(APPEND sources
  ${Pstreams}/Pstream.C
  ${Pstreams}/PstreamCommsStruct.C
  ${Pstreams}/PstreamGlobals.C
  ${Pstreams}/IPstream.C
  ${Pstreams}/OPstream.C
  ${Pstreams}/IPread.C
  ${Pstreams}/OPwrite.C
  ${Pstreams}/PstreamsPrint.C
)

set(dictionary db/dictionary)
list(APPEND sources
  ${dictionary}/dictionary.C
  ${dictionary}/dictionaryIO.C
)

set(entry ${dictionary}/entry)
list(APPEND sources
  ${entry}/entry.C
  ${entry}/entryIO.C
)

set(primitiveEntry ${dictionary}/primitiveEntry)
list(APPEND sources
  ${primitiveEntry}/primitiveEntry.C
  ${primitiveEntry}/primitiveEntryIO.C
)

set(dictionaryEntry ${dictionary}/dictionaryEntry)
list(APPEND sources
  ${dictionaryEntry}/dictionaryEntry.C
  ${dictionaryEntry}/dictionaryEntryIO.C
)

set(functionEntries ${dictionary}/functionEntries)
list(APPEND sources
  ${functionEntries}/functionEntry/functionEntry.C
  ${functionEntries}/includeEntry/includeEntry.C
  ${functionEntries}/includeIfPresentEntry/includeIfPresentEntry.C
  ${functionEntries}/inputModeEntry/inputModeEntry.C
  ${functionEntries}/removeEntry/removeEntry.C
)

set(IOdictionary db/IOobjects/IOdictionary)
list(APPEND sources
  ${IOdictionary}/IOdictionary.C
  ${IOdictionary}/IOdictionaryIO.C
  db/IOobjects/IOMap/IOMapName.C
)

set(IOobject db/IOobject)
list(APPEND sources
  ${IOobject}/IOobject.C
  ${IOobject}/IOobjectIO.C
  ${IOobject}/IOobjectReadHeader.C
  ${IOobject}/IOobjectWriteHeader.C
)

set(regIOobject db/regIOobject)
list(APPEND sources
  ${regIOobject}/regIOobject.C
  ${regIOobject}/regIOobjectRead.C
  ${regIOobject}/regIOobjectWrite.C
  db/IOobjectList/IOobjectList.C
  db/objectRegistry/objectRegistry.C
  db/postfixedSubRegistry/postfixedSubRegistry.C
  db/CallbackRegistry/CallbackRegistryName.C
  db/dlLibraryTable/dlLibraryTable.C
  db/functionObjects/functionObject/functionObject.C
  db/functionObjects/functionObjectList/functionObjectList.C
  db/functionObjects/outputFilterOutputControl/outputFilterOutputControl.C
)

set(Time db/Time)
list(APPEND sources
  ${Time}/TimePaths.C
  ${Time}/TimeState.C
  ${Time}/foamTime.C
  ${Time}/foamTimeIO.C
  ${Time}/findTimes.C
  ${Time}/subCycleTime.C
  ${Time}/findInstance.C
  ${Time}/timeSelector.C
  ${Time}/instant/instant.C
  db/scalarRange/scalarRange.C
  db/scalarRange/scalarRanges.C
  dimensionSet/dimensionSet.C
  dimensionSet/dimensionSetIO.C
  dimensionSet/dimensionSets.C
  dimensionedTypes/dimensionedScalar/dimensionedScalar.C
  dimensionedTypes/dimensionedSphericalTensor/dimensionedSphericalTensor.C
  dimensionedTypes/dimensionedDiagTensor/dimensionedDiagTensor.C
  dimensionedTypes/dimensionedSymmTensor/dimensionedSymmTensor.C
  dimensionedTypes/dimensionedSymmTensor4thOrder/dimensionedSymmTensor4thOrder.C
  dimensionedTypes/dimensionedTensor/dimensionedTensor.C
  dimensionedTypes/dimensionedVectorTensorN/dimensionedVectorTensorN.C
  matrices/solution/solution.C
  matrices/constraint/scalarConstraint.C
)

set(scalarMatrices matrices/scalarMatrices)
list(APPEND sources
  ${scalarMatrices}/scalarMatrices.C
  ${scalarMatrices}/scalarSquareMatrix.C
  ${scalarMatrices}/SVD/SVD.C
)

set(LUscalarMatrix matrices/LUscalarMatrix)
list(APPEND sources
  ${LUscalarMatrix}/LUscalarMatrix.C
  ${LUscalarMatrix}/procLduMatrix.C
  ${LUscalarMatrix}/procLduInterface.C
)

set(lduMatrix matrices/lduMatrix)
list(APPEND sources
  ${lduMatrix}/lduMatrix/lduMatrix.C
  ${lduMatrix}/lduMatrix/lduMatrixOperations.C
  ${lduMatrix}/lduMatrix/lduMatrixATmul.C
  ${lduMatrix}/lduMatrix/lduMatrixUpdateMatrixInterfaces.C
  ${lduMatrix}/lduMatrix/lduMatrixSolver.C
  ${lduMatrix}/lduMatrix/lduMatrixSmoother.C
  ${lduMatrix}/lduMatrix/lduMatrixPreconditioner.C
  ${lduMatrix}/lduMatrix/extendedLduMatrix/extendedLduMatrix.C
  ${lduMatrix}/solvers/diagonalSolver/diagonalSolver.C
  ${lduMatrix}/solvers/smoothSolver/smoothSolver.C
  ${lduMatrix}/solvers/PCG/PCG.C
  ${lduMatrix}/solvers/PBiCG/PBiCG.C
  ${lduMatrix}/solvers/ICCG/ICCG.C
  ${lduMatrix}/solvers/BICCG/BICCG.C
  ${lduMatrix}/smoothers/GaussSeidel/GaussSeidelSmoother.C
  ${lduMatrix}/smoothers/DIC/DICSmoother.C
  ${lduMatrix}/smoothers/DICGaussSeidel/DICGaussSeidelSmoother.C
  ${lduMatrix}/smoothers/DILU/DILUSmoother.C
  ${lduMatrix}/smoothers/DILUGaussSeidel/DILUGaussSeidelSmoother.C
  ${lduMatrix}/preconditioners/noPreconditioner/noPreconditioner.C
  ${lduMatrix}/preconditioners/diagonalPreconditioner/diagonalPreconditioner.C
  ${lduMatrix}/preconditioners/DICPreconditioner/DICPreconditioner.C
  ${lduMatrix}/preconditioners/FDICPreconditioner/FDICPreconditioner.C
  ${lduMatrix}/preconditioners/DILUPreconditioner/DILUPreconditioner.C
  ${lduMatrix}/preconditioners/GAMGPreconditioner/GAMGPreconditioner.C
)

set(lduAddressing ${lduMatrix}/lduAddressing)
list(APPEND sources
  ${lduAddressing}/lduAddressing.C
  ${lduAddressing}/extendedLduAddressing/extendedLduAddressing.C
)

set(lduInterfaces ${lduAddressing}/lduInterfaces)
list(APPEND sources
  ${lduInterfaces}/lduInterface/lduInterface.C
  ${lduInterfaces}/processorLduInterface/processorLduInterface.C
  ${lduInterfaces}/cyclicLduInterface/cyclicLduInterface.C
  ${lduInterfaces}/ggiLduInterface/ggiLduInterface.C
  ${lduInterfaces}/overlapGGILduInterface/overlapGGILduInterface.C
  ${lduInterfaces}/mixingPlaneLduInterface/mixingPlaneLduInterface.C
  ${lduInterfaces}/regionCoupleLduInterface/regionCoupleLduInterface.C
)

set(lduInterfaceFields ${lduAddressing}/lduInterfaceFields)
list(APPEND sources
  ${lduInterfaceFields}/lduInterfaceField/lduInterfaceField.C
  ${lduInterfaceFields}/processorLduInterfaceField/processorLduInterfaceField.C
  ${lduInterfaceFields}/cyclicLduInterfaceField/cyclicLduInterfaceField.C
  ${lduInterfaceFields}/ggiLduInterfaceField/ggiLduInterfaceField.C
  ${lduInterfaceFields}/mixingPlaneLduInterfaceField/mixingPlaneLduInterfaceField.C
  ${lduInterfaceFields}/overlapGGILduInterfaceField/overlapGGILduInterfaceField.C
  ${lduInterfaceFields}/regionCoupleLduInterfaceField/regionCoupleLduInterfaceField.C
)

set(AMG ${lduMatrix}/solvers/AMG)
list(APPEND sources
  ${AMG}/GAMGSolver.C
  ${AMG}/GAMGSolverAgglomerateMatrix.C
  ${AMG}/GAMGSolverScalingFactor.C
  ${AMG}/GAMGSolverSolve.C
)

set(AMGInterfaces ${AMG}/interfaces)
list(APPEND sources
  ${AMGInterfaces}/AMGInterface/AMGInterface.C
  ${AMGInterfaces}/AMGInterface/newAMGInterface.C
  ${AMGInterfaces}/processorAMGInterface/processorAMGInterface.C
  ${AMGInterfaces}/cyclicAMGInterface/cyclicAMGInterface.C
  ${AMGInterfaces}/ggiAMGInterface/ggiAMGInterface.C
  ${AMGInterfaces}/cyclicGGIAMGInterface/cyclicGGIAMGInterface.C
  ${AMGInterfaces}/regionCoupleAMGInterface/regionCoupleAMGInterface.C
  ${AMGInterfaces}/mixingPlaneAMGInterface/mixingPlaneAMGInterface.C
)

set(AMGInterfaceFields ${AMG}/interfaceFields)
list(APPEND sources
  ${AMGInterfaceFields}/AMGInterfaceField/AMGInterfaceField.C
  ${AMGInterfaceFields}/AMGInterfaceField/newAMGInterfaceField.C
  ${AMGInterfaceFields}/processorAMGInterfaceField/processorAMGInterfaceField.C
  ${AMGInterfaceFields}/cyclicAMGInterfaceField/cyclicAMGInterfaceField.C
  ${AMGInterfaceFields}/ggiAMGInterfaceField/ggiAMGInterfaceField.C
  ${AMGInterfaceFields}/cyclicGgiAMGInterfaceField/cyclicGgiAMGInterfaceField.C
  ${AMGInterfaceFields}/regionCoupleAMGInterfaceField/regionCoupleAMGInterfaceField.C
  ${AMGInterfaceFields}/mixingPlaneAMGInterfaceField/mixingPlaneAMGInterfaceField.C
)

set(AMGAgglomerations ${AMG}/AMGAgglomerations)

set(GAMGAgglomeration ${AMGAgglomerations}/GAMGAgglomeration)
list(APPEND sources
  ${GAMGAgglomeration}/GAMGAgglomeration.C
  ${GAMGAgglomeration}/GAMGAgglomerateLduAddressing.C
)

set(pairGAMGAgglomeration ${AMGAgglomerations}/pairGAMGAgglomeration)
list(APPEND sources
  ${pairGAMGAgglomeration}/pairGAMGAgglomeration.C
  ${pairGAMGAgglomeration}/pairGAMGAgglomerate.C
  ${pairGAMGAgglomeration}/pairGAMGAgglomerationCombineLevels.C
)

set(algebraicPairGAMGAgglomeration ${AMGAgglomerations}/algebraicPairGAMGAgglomeration)
list(APPEND sources
  ${algebraicPairGAMGAgglomeration}/algebraicPairGAMGAgglomeration.C
  meshes/lduMesh/lduMesh.C
  meshes/MeshObject/meshObjectBase.C
)

set(primitiveShapes meshes/primitiveShapes)

list(APPEND sources
  ${primitiveShapes}/line/line.C
  ${primitiveShapes}/plane/plane.C
  ${primitiveShapes}/triangle/triangleFuncs.C
  ${primitiveShapes}/triangle/intersection.C
)

set(meshShapes meshes/meshShapes)
list(APPEND sources
  ${meshShapes}/edge/edge.C
  ${meshShapes}/edge/edgeIOList.C
)

set(face ${meshShapes}/face)
list(APPEND sources
  ${face}/face.C
  ${face}/faceIntersection.C
  ${face}/faceContactSphere.C
  ${face}/faceAreaInContact.C
  ${face}/faceIOList.C
)

set(cell ${meshShapes}/cell)
list(APPEND sources
  ${cell}/cell.C
  ${cell}/oppositeCellFace.C
  ${cell}/cellIOList.C
)

set(tetCell ${meshShapes}/tetCell)
list(APPEND sources
  ${tetCell}/tetCell.C
)

set(cellModeller ${meshShapes}/cellModeller)
list(APPEND sources
  ${cellModeller}/cellModeller.C
)

set(cellModel ${meshShapes}/cellModel)
list(APPEND sources
  ${cellModel}/cellModel.C
  ${cellModel}/cellModelIO.C
)

set(cellShape ${meshShapes}/cellShape)
list(APPEND sources
  ${cellShape}/cellShape.C
  ${cellShape}/cellShapeEqual.C
  ${cellShape}/cellShapeIO.C
  ${cellShape}/cellShapeIOList.C
  meshes/patchIdentifier/patchIdentifier.C
)

set(polyMesh meshes/polyMesh)

set(polyPatches ${polyMesh}/polyPatches)
list(APPEND sources
  ${polyPatches}/polyPatch/polyPatch.C
  ${polyPatches}/polyPatch/newPolyPatch.C
)

set(basicPolyPatches ${polyPatches}/basic)
list(APPEND sources
  ${basicPolyPatches}/coupled/coupledPolyPatch.C
  ${basicPolyPatches}/generic/genericPolyPatch.C
)

set(constraintPolyPatches ${polyPatches}/constraint)
list(APPEND sources
  ${constraintPolyPatches}/empty/emptyPolyPatch.C
  ${constraintPolyPatches}/symmetry/symmetryPolyPatch.C
  ${constraintPolyPatches}/wedge/wedgePolyPatch.C
  ${constraintPolyPatches}/cyclic/cyclicPolyPatch.C
  ${constraintPolyPatches}/processor/processorPolyPatch.C
  ${constraintPolyPatches}/ggi/ggiPolyPatch.C
  ${constraintPolyPatches}/cyclicGgi/cyclicGgiPolyPatch.C
  ${constraintPolyPatches}/overlapGgi/overlapGgiPolyPatch.C
  ${constraintPolyPatches}/overlapGgi/overlapGgiPolyPatchGeometry.C
  ${constraintPolyPatches}/mixingPlane/mixingPlanePolyPatch.C
  ${constraintPolyPatches}/regionCouple/regionCouplePolyPatch.C
)

set(derivedPolyPatches ${polyPatches}/derived)
list(APPEND sources
  ${derivedPolyPatches}/wall/wallPolyPatch.C
  ${derivedPolyPatches}/cohesive/cohesivePolyPatch.C
)

set(polyBoundaryMesh ${polyMesh}/polyBoundaryMesh)
list(APPEND sources
  ${polyBoundaryMesh}/polyBoundaryMesh.C
  ${polyBoundaryMesh}/polyBoundaryMeshEntries.C
  meshes/ProcessorTopology/commSchedule.C
)

set(globalMeshData ${polyMesh}/globalMeshData)
list(APPEND sources
  ${globalMeshData}/globalMeshData.C
  ${globalMeshData}/globalPoints.C
  ${globalMeshData}/globalIndex.C
  ${polyMesh}/syncTools/syncTools.C
)

set(zones ${polyMesh}/zones)

set(cellZone ${zones}/cellZone)
list(APPEND sources
  ${cellZone}/cellZone.C
  ${cellZone}/newCellZone.C
)

set(faceZone ${zones}/faceZone)
list(APPEND sources
  ${faceZone}/faceZone.C
  ${faceZone}/newFaceZone.C
)

set(pointZone ${polyMesh}/zones/pointZone)
list(APPEND sources
  ${pointZone}/pointZone.C
  ${pointZone}/newPointZone.C
  ${polyMesh}/polyMesh.C
  ${polyMesh}/polyMeshFromShapeMesh.C
  ${polyMesh}/polyMeshIO.C
  ${polyMesh}/polyMeshInitMesh.C
  ${polyMesh}/polyMeshClear.C
  ${polyMesh}/polyMeshUpdate.C
)

set(primitiveMesh meshes/primitiveMesh)
list(APPEND sources
  ${primitiveMesh}/primitiveMesh.C
  ${primitiveMesh}/primitiveMeshCellCells.C
  ${primitiveMesh}/primitiveMeshCellCentresAndVols.C
  ${primitiveMesh}/primitiveMeshCellEdges.C
  ${primitiveMesh}/primitiveMeshCells.C
  ${primitiveMesh}/primitiveMeshClear.C
  ${primitiveMesh}/primitiveMeshEdgeCells.C
  ${primitiveMesh}/primitiveMeshEdgeFaces.C
  ${primitiveMesh}/primitiveMeshEdges.C
  ${primitiveMesh}/primitiveMeshFaceCentresAndAreas.C
  ${primitiveMesh}/primitiveMeshFindCell.C
  ${primitiveMesh}/primitiveMeshPointCells.C
  ${primitiveMesh}/primitiveMeshPointFaces.C
  ${primitiveMesh}/primitiveMeshPointEdges.C
  ${primitiveMesh}/primitiveMeshPointPoints.C
  ${primitiveMesh}/primitiveMeshCellPoints.C
  ${primitiveMesh}/primitiveMeshCalcCellShapes.C
)

set(primitiveMeshCheck ${primitiveMesh}/primitiveMeshCheck)
list(APPEND sources
  ${primitiveMeshCheck}/primitiveMeshCheck.C
  ${primitiveMeshCheck}/primitiveMeshCheckMotion.C
  ${primitiveMeshCheck}/primitiveMeshCheckPointNearness.C
  ${primitiveMeshCheck}/primitiveMeshCheckEdgeLength.C
)

set(primitivePatch ${primitiveMesh}/primitivePatch)
list(APPEND sources
  ${primitivePatch}/patchZones.C
  ${primitivePatch}/walkPatch.C
)

set(cellMatcher ${meshShapes}/cellMatcher)
list(APPEND sources
  ${cellMatcher}/cellMatcher.C
  ${cellMatcher}/hexMatcher.C
  ${cellMatcher}/wedgeMatcher.C
  ${cellMatcher}/prismMatcher.C
  ${cellMatcher}/pyrMatcher.C
  ${cellMatcher}/tetWedgeMatcher.C
  ${cellMatcher}/tetMatcher.C
  ${cellMatcher}/degenerateMatcher.C
)

set(mapPolyMesh ${polyMesh}/mapPolyMesh)
list(APPEND sources
  ${mapPolyMesh}/mapPolyMesh.C
  ${mapPolyMesh}/pointMapper/pointMapper.C
  ${mapPolyMesh}/faceMapper/faceMapper.C
  ${mapPolyMesh}/cellMapper/cellMapper.C
  ${mapPolyMesh}/mapDistribute/mapDistribute.C
  ${mapPolyMesh}/mapDistribute/mapDistributePolyMesh.C
  ${mapPolyMesh}/mapAddedPolyMesh.C
)

set(PrimitivePatch ${primitiveMesh}/PrimitivePatchTemplate)
list(APPEND sources
  ${PrimitivePatch}/PrimitivePatchName.C
)

set(pointMesh meshes/pointMesh)
list(APPEND sources
  ${pointMesh}/pointMesh.C
)

set(pointMeshMapper ${pointMesh}/pointMeshMapper)
list(APPEND sources
  ${pointMeshMapper}/pointPatchMapper.C
)

set(pointPatches ${pointMesh}/pointPatches)
list(APPEND sources
  ${pointPatches}/pointPatch/pointPatch.C
  ${pointPatches}/facePointPatch/facePointPatch.C
  ${pointPatches}/facePointPatch/newFacePointPatch.C
)

set(basicPointPatches ${pointPatches}/basic)
list(APPEND sources
  ${basicPointPatches}/coupled/coupledPointPatch.C
  ${basicPointPatches}/generic/genericPointPatch.C
)

set(constraintPointPatches ${pointPatches}/constraint)
list(APPEND sources
  ${constraintPointPatches}/empty/emptyPointPatch.C
  ${constraintPointPatches}/symmetry/symmetryPointPatch.C
  ${constraintPointPatches}/wedge/wedgePointPatch.C
  ${constraintPointPatches}/cyclic/cyclicPointPatch.C
  ${constraintPointPatches}/processor/processorPointPatch.C
  ${constraintPointPatches}/ggi/ggiPointPatch.C
  ${constraintPointPatches}/cyclicGgi/cyclicGgiPointPatch.C
  ${constraintPointPatches}/overlapGgi/overlapGgiPointPatch.C
  ${constraintPointPatches}/mixingPlane/mixingPlanePointPatch.C
  ${constraintPointPatches}/regionCouple/regionCouplePointPatch.C
)

set(derivedPointPatches ${pointPatches}/derived)
list(APPEND sources
  ${derivedPointPatches}/coupled/coupledFacePointPatch.C
  ${derivedPointPatches}/global/globalPointPatch.C
  ${derivedPointPatches}/wall/wallPointPatch.C
  ${derivedPointPatches}/cohesive/cohesivePointPatch.C
)

set(pointBoundaryMesh ${pointMesh}/pointBoundaryMesh)
list(APPEND sources
  ${pointBoundaryMesh}/pointBoundaryMesh.C
  meshes/boundBox/boundBox.C
)

set(meshTools meshes/meshTools)
list(APPEND sources
  ${meshTools}/meshTools.C
  ${meshTools}/matchPoints.C
  ${meshTools}/mergePoints.C
  fields/UniformDimensionedFields/uniformDimensionedFields.C
  fields/cloud/cloud.C
)

set(Fields fields/Fields)
list(APPEND sources
  ${Fields}/labelField/labelField.C
  ${Fields}/scalarField/scalarField.C
  ${Fields}/vectorField/vectorField.C
  ${Fields}/sphericalTensorField/sphericalTensorField.C
  ${Fields}/diagTensorField/diagTensorField.C
  ${Fields}/symmTensorField/symmTensorField.C
  ${Fields}/symmTensor4thOrderField/symmTensor4thOrderField.C
  ${Fields}/tensorField/tensorField.C
  ${Fields}/complexFields/complexFields.C
  ${Fields}/VectorNFields/VectorNFields.C
  ${Fields}/VectorNFields/SphericalTensorNFields.C
  ${Fields}/VectorNFields/DiagTensorNFields.C
  ${Fields}/VectorNFields/TensorNFields.C
  ${Fields}/VectorNFields/expandContract/ExpandTensorNField.C
  ${Fields}/labelField/labelIOField.C
  ${Fields}/scalarField/scalarIOField.C
  ${Fields}/vectorField/vectorIOField.C
  ${Fields}/vector2DField/vector2DIOField.C
  ${Fields}/sphericalTensorField/sphericalTensorIOField.C
  ${Fields}/diagTensorField/diagTensorIOField.C
  ${Fields}/symmTensorField/symmTensorIOField.C
  ${Fields}/symmTensor4thOrderField/symmTensor4thOrderIOField.C
  ${Fields}/tensorField/tensorIOField.C
  ${Fields}/transformField/transformField.C
)

set(pointPatchFields fields/pointPatchFields)
list(APPEND sources
  ${pointPatchFields}/pointPatchField/pointPatchFields.C
)

set(basicPointPatchFields ${pointPatchFields}/basic)
list(APPEND sources
  ${basicPointPatchFields}/calculated/calculatedPointPatchFields.C
  ${basicPointPatchFields}/generic/genericPointPatchFields.C
  ${basicPointPatchFields}/coupled/coupledPointPatchFields.C
  ${basicPointPatchFields}/value/valuePointPatchFields.C
  ${basicPointPatchFields}/fixedValue/fixedValuePointPatchFields.C
  ${basicPointPatchFields}/zeroGradient/zeroGradientPointPatchFields.C
  ${basicPointPatchFields}/mixed/mixedPointPatchFields.C
)

set(constraintPointPatchFields ${pointPatchFields}/constraint)
list(APPEND sources
  ${constraintPointPatchFields}/empty/emptyPointPatchFields.C
  ${constraintPointPatchFields}/symmetry/symmetryPointPatchFields.C
  ${constraintPointPatchFields}/wedge/wedgePointPatchFields.C
  ${constraintPointPatchFields}/cyclic/cyclicPointPatchFields.C
  ${constraintPointPatchFields}/processor/processorPointPatchFields.C
)

set(derivedPointPatchFields ${pointPatchFields}/derived)
list(APPEND sources
  ${derivedPointPatchFields}/slip/slipPointPatchFields.C
  ${derivedPointPatchFields}/global/globalPointPatchFields.C
  ${derivedPointPatchFields}/uniformFixedValue/uniformFixedValuePointPatchFields.C
  ${derivedPointPatchFields}/timeVaryingUniformFixedValue/timeVaryingUniformFixedValuePointPatchFields.C
  ${derivedPointPatchFields}/oscillatingFixedValue/oscillatingFixedValuePointPatchFields.C
  fields/GeometricFields/pointFields/pointFields.C
  meshes/bandCompression/bandCompression.C
  meshes/preservePatchTypes/preservePatchTypes.C
)

set(interpolations interpolations)
set(interpolation ${interpolations}/interpolation)
list(APPEND sources
  ${interpolations}/splineInterpolateXY/splineInterpolateXY.C
  ${interpolations}/patchToPatchInterpolation/PatchToPatchInterpolationName.C
  ${interpolations}/GGIInterpolation/GGIInterpolationName.C
  ${interpolations}/MixingPlaneInterpolation/MixingPlaneInterpolationName.C
  ${interpolations}/RBFInterpolation/RBFInterpolation.C
  ${interpolations}/RBFInterpolation/RBFFunctions/RBFFunction/RBFFunction.C
  ${interpolations}/RBFInterpolation/RBFFunctions/RBFFunction/newRBFFunction.C
  ${interpolations}/RBFInterpolation/RBFFunctions/W2/W2.C
  ${interpolations}/RBFInterpolation/RBFFunctions/Gauss/Gauss.C
  ${interpolations}/RBFInterpolation/RBFFunctions/TPS/TPS.C
  ${interpolations}/RBFInterpolation/RBFFunctions/IMQB/IMQB.C
  algorithms/MeshWave/MeshWaveName.C
  algorithms/MeshWave/FaceCellWaveName.C
  algorithms/polygon/clipping/SutherlandHodgman.C
  algorithms/polygon/pointInPolygon/HormannAgathos.C
  algorithms/rotation/RodriguesRotation.C
)

set(octree algorithms/octree/octree)
list(APPEND sources
  ${octree}/octreeName.C
  ${octree}/octreeDataPoint.C
  ${octree}/octreeDataPointTreeLeaf.C
  ${octree}/octreeDataEdges.C
  ${octree}/octreeDataCell.C
  ${octree}/octreeDataFace.C
  ${octree}/octreeDataBoundBox.C
  ${octree}/treeBoundBox.C
  ${octree}/treeNodeName.C
  ${octree}/treeLeafName.C
  ${octree}/pointIndexHitIOList.C
)

set(indexedOctree algorithms/octree/indexedOctree)
list(APPEND sources
  ${indexedOctree}/indexedOctreeName.C
  ${indexedOctree}/treeDataCell.C
  ${indexedOctree}/treeDataEdge.C
  ${indexedOctree}/treeDataFace.C
  ${indexedOctree}/treeDataPoint.C
  graph/curve/curve.C
  graph/graph.C
)

set(writers graph/writers)
list(APPEND sources
  ${writers}/rawGraph/rawGraph.C
  ${writers}/gnuplotGraph/gnuplotGraph.C
  ${writers}/xmgrGraph/xmgrGraph.C
  ${writers}/jplotGraph/jplotGraph.C
  primitives/BlockCoeff/blockCoeffBase.C
  primitives/BlockCoeff/scalarBlockCoeff.C
  primitives/BlockCoeff/sphericalTensorBlockCoeff.C
  primitives/BlockCoeff/symmTensorBlockCoeff.C
  primitives/BlockCoeff/symmTensor4thOrderBlockCoeff.C
  primitives/BlockCoeff/tensorBlockCoeff.C
)

set(BlockCoeffNorm primitives/BlockCoeff/BlockCoeffNorm)
list(APPEND sources
  ${BlockCoeffNorm}/BlockCoeffNorm/blockCoeffNorms.C
  ${BlockCoeffNorm}/BlockCoeffTwoNorm/blockCoeffTwoNorms.C
  ${BlockCoeffNorm}/BlockCoeffMaxNorm/blockCoeffMaxNorms.C
  ${BlockCoeffNorm}/BlockCoeffComponentNorm/blockCoeffComponentNorms.C
  fields/expandContract/expandTensorField.C
  fields/CoeffField/scalarCoeffField.C
  fields/CoeffField/sphericalTensorCoeffField.C
  fields/CoeffField/diagTensorCoeffField.C
  fields/CoeffField/symmTensorCoeffField.C
  fields/CoeffField/symmTensor4thOrderCoeffField.C
  fields/CoeffField/tensorCoeffField.C
)

set(BlockLduMatrix matrices/blockLduMatrix/BlockLduMatrix)
list(APPEND sources
  ${BlockLduMatrix}/blockLduMatrixName.C
  ${BlockLduMatrix}/scalarBlockLduMatrix.C
  ${BlockLduMatrix}/sphericalTensorBlockLduMatrix.C
  ${BlockLduMatrix}/symmTensorBlockLduMatrix.C
  ${BlockLduMatrix}/tensorBlockLduMatrix.C
  ${BlockLduMatrix}/BlockConstraint/scalarBlockConstraint.C
  ${BlockLduMatrix}/extendedBlockLduMatrix/extendedBlockLduMatrices.C
  ${BlockLduMatrix}/extendedBlockLduMatrix/sphericalTensorExtendedBlockLduMatrix.C
  ${BlockLduMatrix}/extendedBlockLduMatrix/symmTensorExtendedBlockLduMatrix.C
  ${BlockLduMatrix}/extendedBlockLduMatrix/tensorExtendedBlockLduMatrix.C
  ${BlockLduMatrix}/extendedBlockLduMatrix/extendedBlockVectorNMatrices.C
)

set(BlockLduInterfaceFields ${BlockLduMatrix}/BlockLduInterfaceFields)
list(APPEND sources
  ${BlockLduInterfaceFields}/BlockLduInterfaceField/BlockLduInterfaceFields.C
  ${BlockLduInterfaceFields}/BlockLduInterfaceField/blockVectorNLduInterfaceFields.C
  ${BlockLduInterfaceFields}/GGIBlockLduInterfaceField/GGIBlockLduInterfaceFields.C
  ${BlockLduInterfaceFields}/MixingPlaneBlockLduInterfaceField/MixingPlaneBlockLduInterfaceFields.C
  ${BlockLduInterfaceFields}/OverlapGGIBlockLduInterfaceField/OverlapGGIBlockLduInterfaceFields.C
  ${BlockLduInterfaceFields}/ProcessorBlockLduInterfaceField/ProcessorBlockLduInterfaceFields.C
)

set(BlockAMG matrices/blockLduMatrix/BlockAMG)
list(APPEND sources
  ${BlockAMG}/blockAMGCycles.C
  ${BlockAMG}/blockAMGLevels.C
)

set(BlockAMGInterfaceFields ${BlockAMG}/BlockAMGInterfaceFields)
list(APPEND sources
  ${BlockAMGInterfaceFields}/BlockAMGInterfaceField/blockAMGInterfaceFields.C
  ${BlockAMGInterfaceFields}/BlockAMGInterfaceField/blockVectorNAMGInterfaceFields.C
  ${BlockAMGInterfaceFields}/ProcessorBlockAMGInterfaceField/ProcessorBlockAMGInterfaceFields.C
  ${BlockAMGInterfaceFields}/GGIBlockAMGInterfaceField/GGIBlockAMGInterfaceFields.C
)

set(BlockMatrixCoarsening ${BlockAMG}/BlockMatrixCoarsening)
list(APPEND sources
  ${BlockMatrixCoarsening}/BlockMatrixCoarsening/blockMatrixCoarsenings.C
  ${BlockMatrixCoarsening}/BlockMatrixAgglomeration/blockMatrixAgglomerations.C
)

set(BlockLduPrecons matrices/blockLduMatrix/BlockLduPrecons)
list(APPEND sources
  ${BlockLduPrecons}/BlockLduPrecon/blockLduPrecons.C
  ${BlockLduPrecons}/BlockNoPrecon/blockNoPrecons.C
  ${BlockLduPrecons}/BlockDiagonalPrecon/scalarBlockDiagonalPrecon.C
  ${BlockLduPrecons}/BlockDiagonalPrecon/tensorBlockDiagonalPrecon.C
  ${BlockLduPrecons}/BlockDiagonalPrecon/blockDiagonalPrecons.C
  ${BlockLduPrecons}/BlockGaussSeidelPrecon/scalarBlockGaussSeidelPrecon.C
  ${BlockLduPrecons}/BlockGaussSeidelPrecon/tensorBlockGaussSeidelPrecon.C
  ${BlockLduPrecons}/BlockGaussSeidelPrecon/blockGaussSeidelPrecons.C
  ${BlockLduPrecons}/BlockCholeskyPrecon/scalarBlockCholeskyPrecon.C
  ${BlockLduPrecons}/BlockCholeskyPrecon/tensorBlockCholeskyPrecon.C
  ${BlockLduPrecons}/BlockCholeskyPrecon/blockCholeskyPrecons.C
  ${BlockLduPrecons}/BlockAMGPrecon/blockAMGPrecons.C
  ${BlockLduPrecons}/BlockILUCpPrecon/scalarBlockILUCpPrecon.C
  ${BlockLduPrecons}/BlockILUCpPrecon/tensorBlockILUCpPrecon.C
  ${BlockLduPrecons}/BlockILUCpPrecon/blockILUCpPrecons.C
)

set(BlockLduSmoothers matrices/blockLduMatrix/BlockLduSmoothers)
list(APPEND sources
  ${BlockLduSmoothers}/BlockLduSmoother/blockLduSmoothers.C
  ${BlockLduSmoothers}/BlockGaussSeidelSmoother/blockGaussSeidelSmoothers.C
  ${BlockLduSmoothers}/BlockILUSmoother/blockILUSmoothers.C
  ${BlockLduSmoothers}/BlockILUCpSmoother/blockILUCpSmoothers.C
)

set(BlockLduSolvers matrices/blockLduMatrix/BlockLduSolvers)
list(APPEND sources
  ${BlockLduSolvers}/blockVectorNSolvers.C
  ${BlockLduSolvers}/BlockLduSolver/blockLduSolvers.C
  ${BlockLduSolvers}/BlockDiagonal/blockDiagonalSolvers.C
  ${BlockLduSolvers}/BlockGaussSeidel/blockGaussSeidelSolvers.C
  ${BlockLduSolvers}/BlockCG/blockCGSolvers.C
  ${BlockLduSolvers}/BlockBiCGStab/blockBiCGStabSolvers.C
  ${BlockLduSolvers}/BlockGMRES/blockGMRESSolvers.C
  ${BlockLduSolvers}/BlockAMGSolver/blockAMGSolvers.C
  ${BlockLduSolvers}/Segregated/segregatedSolvers.C
)

add_foam_library(foam SHARED ${sources})

# This is a prototype! Replace second argument with a list of the highest
# dependencies only. Inherited dependencies will be treated automatically.
target_link_libraries(foam OSspecific mpi ZLIB::ZLIB)
add_dependencies(foam getGitVersion)
