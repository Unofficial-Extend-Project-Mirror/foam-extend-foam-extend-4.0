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

set(meshSurfaceEngine utilities/surfaceTools/meshSurfaceEngine)
set(meshSurfacePartitioner utilities/surfaceTools/meshSurfacePartitioner)

set(bndLayers utilities/boundaryLayers)
set(boundaryLayers ${bndLayers}/boundaryLayers)
set(detectBoundaryLayers ${bndLayers}/detectBoundaryLayers)
set(extrudeLayer ${bndLayers}/extrudeLayer)
set(refineBoundaryLayers ${bndLayers}/refineBoundaryLayers)
set(triangulateNonPlanarBaseFaces ${bndLayers}/triangulateNonPlanarBaseFaces)

set(checkMeshDict utilities/checkMeshDict)

set(anisotropicMeshing utilities/anisotropicMeshing)
set(coordinateModification ${anisotropicMeshing}/coordinateModification)
set(edgeMeshGeometryModification ${anisotropicMeshing}/edgeMeshGeometryModification)
set(surfaceMeshGeometryModification ${anisotropicMeshing}/surfaceMeshGeometryModification)
set(polyMeshGenGeometryModification ${anisotropicMeshing}/polyMeshGenGeometryModification)

set(meshSurfaceCheckInvertedVertices utilities/surfaceTools/meshSurfaceCheckInvertedVertices)
set(meshSurfaceCheckEdgeTypes utilities/surfaceTools/meshSurfaceCheckEdgeTypes)
set(meshSurfaceCutter utilities/surfaceTools/meshSurfaceCutter)
set(meshSurfaceMapper utilities/surfaceTools/meshSurfaceMapper)
set(meshSurfaceMapper2D utilities/surfaceTools/meshSurfaceMapper2D)
set(edgeExtraction utilities/surfaceTools/edgeExtraction)
set(edgeExtractor ${edgeExtraction}/edgeExtractor)
set(meshSurfaceEdgeExtractor utilities/surfaceTools/meshSurfaceEdgeExtractor)
set(meshSurfaceEdgeExtractorNonTopo utilities/surfaceTools/meshSurfaceEdgeExtractorNonTopo)
set(meshSurfaceEdgeExtractor2D utilities/surfaceTools/meshSurfaceEdgeExtractor2D)
set(meshSurfaceEdgeExtractorFUN utilities/surfaceTools/meshSurfaceEdgeExtractorFUN)

set(createFundamentalSheetsBase utilities/surfaceTools/createFundamentalSheets)
set(createFundamentalSheets ${createFundamentalSheetsBase}/createFundamentalSheets)
set(createFundamentalSheetsFJ ${createFundamentalSheetsBase}/createFundamentalSheetsFJ)
set(createFundamentalSheetsJFS ${createFundamentalSheetsBase}/createFundamentalSheetsJFS)

set(correctEdgesBetweenPatches utilities/surfaceTools/correctEdgesBetweenPatches)

set(decomposeCellsNearConcaveEdges utilities/surfaceTools/decomposeCellsNearConcaveEdges)
set(renameBoundaryPatches utilities/surfaceTools/renameBoundaryPatches)

set(intersectionTools utilities/intersectionTools)
set(findCellsIntersectingSurface ${intersectionTools}/findCellsIntersectingSurface)

set(meshOptimizer utilities/smoothers/geometry/meshOptimizer)
set(tetMeshOptimisation ${meshOptimizer}/tetMeshOptimisation)
set(boundaryLayerOptimisation ${meshOptimizer}/boundaryLayerOptimisation)

set(symmetryPlaneOptimisation ${meshOptimizer}/symmetryPlaneOptimisation)
set(simplexSmoother ${tetMeshOptimisation}/advancedSmoothers/simplexSmoother)
set(knuppMetric ${tetMeshOptimisation}/advancedSmoothers/knuppMetric)
set(meshUntangler ${tetMeshOptimisation}/advancedSmoothers/meshUntangler)
set(quadricMetric ${tetMeshOptimisation}/advancedSmoothers/quadricMetric)
set(volumeOptimizer ${tetMeshOptimisation}/advancedSmoothers/volumeOptimizer)

set(meshSurfaceOptimizer utilities/smoothers/geometry/meshSurfaceOptimizer)
set(surfaceOptimizer ${meshSurfaceOptimizer}/advancedSurfaceSmoothers/surfaceOptimizer)
set(surfaceMorpherCells utilities/surfaceTools/surfaceMorpherCells)

set(utilitiesOctrees utilities/octrees)
set(meshOctree utilities/octrees/meshOctree)
set(meshOctreeCube utilities/octrees/meshOctree/meshOctreeCube)
set(meshOctreeInsideOutside utilities/octrees/meshOctree/meshOctreeInsideOutside)
set(meshOctreeCreator utilities/octrees/meshOctree/meshOctreeCreator)
set(meshOctreeAddressing utilities/octrees/meshOctree/meshOctreeAddressing)
set(meshOctreeModifier utilities/octrees/meshOctree/meshOctreeModifier)
set(meshOctreeAutoRef utilities/octrees/meshOctree/meshOctreeAutomaticRefinement)
set(patchRefinement utilities/octrees/meshOctree/refinementControls/patchRefinement)
set(objectRefinement utilities/octrees/meshOctree/refinementControls/objectRefinement)

set(helperFunctions utilities/helperFunctions)
set(createFacesFromChain utilities/helperClasses/createFacesFromChain)
set(sortEdgesIntoChains utilities/helperClasses/sortEdgesIntoChains)
set(trianglePlaneIntersections utilities/helperClasses/trianglePlaneIntersections)

set(tetCreatorOctree utilities/tetrahedra/tetCreatorOctree)
set(faceDecomposition utilities/faceDecomposition)
set(decomposeCells utilities/decomposeCells)

set(topology utilities/smoothers/topology)
set(topologicalCleaner ${topology}/topologicalCleaner)
set(checkBoundaryFacesSharingTwoEdges ${topology}/checkBoundaryFacesSharingTwoEdges)
set(checkCellConnectionsOverFaces ${topology}/checkCellConnectionsOverFaces)
set(checkIrregularSurfaceConnections ${topology}/checkIrregularSurfaceConnections)
set(checkNonMappableCellConnections ${topology}/checkNonMappableCellConnections)

set(triSurfaceTools utilities/triSurfaceTools)

set(triSurface2DCheck ${triSurfaceTools}/triSurface2DCheck)
set(triSurfaceChecks ${triSurfaceTools}/triSurfaceChecks)
set(triSurfaceCleanupDuplicates ${triSurfaceTools}/triSurfaceCleanupDuplicates)
set(triSurfaceCleanupDuplicateTriangles ${triSurfaceTools}/triSurfaceCleanupDuplicateTriangles)
set(triSurfaceCopyParts ${triSurfaceTools}/triSurfaceCopyParts)
set(triSurfaceCurvatureEstimator ${triSurfaceTools}/triSurfaceCurvatureEstimator)
set(triSurfacePartitioner ${triSurfaceTools}/triSurfacePartitioner)
set(triSurfaceDetectFeatureEdges ${triSurfaceTools}/triSurfaceDetectFeatureEdges)
set(triSurfaceClassifyEdges ${triSurfaceTools}/triSurfaceClassifyEdges)
set(triSurfaceImportSurfaceAsSubset ${triSurfaceTools}/triSurfaceImportSurfaceAsSubset)
set(triSurfacePatchManipulator ${triSurfaceTools}/triSurfacePatchManipulator)
set(triSurfaceRemoveFacets ${triSurfaceTools}/triSurfaceRemoveFacets)
set(triSurfaceExtrude2DEdges ${triSurfaceTools}/triSurfaceExtrude2DEdges)
set(triSurfaceMetaData ${triSurfaceTools}/triSurfaceMetaData)

set(polyMeshGen utilities/meshes/polyMeshGen)
set(boundaryPatch utilities/meshes/polyMeshGen/boundaryPatch)
set(polyMeshGen2DEngine utilities/meshes/polyMeshGen2DEngine)
set(polyMeshGenModifier utilities/meshes/polyMeshGenModifier)
set(polyMeshGenAddressing utilities/meshes/polyMeshGenAddressing)
set(polyMeshGenChecks utilities/meshes/polyMeshGenChecks)
set(partTetMesh utilities/meshes/partTetMesh)
set(partTriMesh utilities/meshes/partTriMesh)
set(primitiveMesh utilities/meshes/primitiveMesh)
set(triSurf utilities/meshes/triSurf)
set(cell utilities/meshes/primitives/cell)
set(edge utilities/meshes/primitives/edge)
set(face utilities/meshes/primitives/face)

set(containers utilities/containers)
set(VRWGraph ${containers}/VRWGraph)
set(VRWGraphList ${containers}/VRWGraphList)
set(graphs ${containers}/Graphs)
set(lists ${containers}/Lists)

set(meshZipper utilities/meshZipper)

set(writeAsFPMA utilities/dataConversion/foamToFPMA)

set(polyMeshExtractor pMeshLibrary/polyMeshExtractor)
set(polyMeshGenerator pMeshLibrary/polyMeshGenerator)

set(cartesianMeshExtractor cartesianMesh/cartesianMeshExtractor)
set(cartesianMeshGenerator cartesianMesh/cartesianMeshGenerator)

set(cartesian2DMeshGenerator cartesian2DMesh/cartesian2DMeshGenerator)

set(tetMeshExtractor tetMesh/tetMeshExtractor)
set(tetMeshExtractorOctree tetMesh/tetMeshExtractorOctree)
set(tetMeshGenerator tetMesh/tetMeshGenerator)

set(voronoiMeshExtractor voronoiMesh/voronoiMeshExtractor)
set(voronoiMeshGenerator voronoiMesh/voronoiMeshGenerator)

set(workflowControls utilities/workflowControls)

list(APPEND SOURCES
  ${checkMeshDict}/checkMeshDict.C
  ${lists}/pointFieldPMG.C
  ${lists}/faceListPMG.C
  ${VRWGraph}/VRWGraph.C
  ${VRWGraph}/VRWGraphSMPModifier.C
  ${VRWGraphList}/VRWGraphList.C
  ${graphs}/cellIOGraph.C
  ${graphs}/faceIOGraph.C
  ${polyMeshGen}/polyMeshGen.C
  ${polyMeshGen}/polyMeshGenPoints.C
  ${polyMeshGen}/polyMeshGenFaces.C
  ${polyMeshGen}/polyMeshGenCells.C
  ${polyMeshGen2DEngine}/polyMeshGen2DEngine.C
  ${coordinateModification}/coordinateModification.C
  ${coordinateModification}/newCoordinateModification.C
  ${coordinateModification}/planeScaling.C
  ${coordinateModification}/boxScaling.C
  ${edgeMeshGeometryModification}/edgeMeshGeometryModification.C
  ${surfaceMeshGeometryModification}/surfaceMeshGeometryModification.C
  ${polyMeshGenGeometryModification}/polyMeshGenGeometryModification.C
  ${coordinateModification}/coordinateModifier.C
  ${boundaryPatch}/boundaryPatchBase.C
  ${boundaryPatch}/boundaryPatch.C
  ${boundaryPatch}/processorBoundaryPatch.C
  ${polyMeshGenModifier}/polyMeshGenModifierRemoveUnusedVertices.C
  ${polyMeshGenModifier}/polyMeshGenModifierRemoveFaces.C
  ${polyMeshGenModifier}/polyMeshGenModifierRemoveCells.C
  ${polyMeshGenModifier}/polyMeshGenModifierReorderBoundaryFaces.C
  ${polyMeshGenModifier}/polyMeshGenModifierAddCells.C
  ${polyMeshGenModifier}/polyMeshGenModifierAddProcessorFaces.C
  ${polyMeshGenModifier}/polyMeshGenModifierAddBufferCells.C
  ${polyMeshGenModifier}/polyMeshGenModifierReplaceBoundary.C
  ${polyMeshGenModifier}/polyMeshGenModifierZipUpCells.C
  ${polyMeshGenModifier}/polyMeshGenModifierRenumberMesh.C
  ${polyMeshGenModifier}/polyMeshGenModifierAddCellByCell.C
  ${polyMeshGenAddressing}/polyMeshGenAddressing.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingCellCells.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingCellEdges.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingCellPoints.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingCentresAndAreas.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingCentresAndVols.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingClear.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingUpdateGeometry.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingEdgeCells.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingEdgeFaces.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingFaceEdges.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingEdges.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingPointCells.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingPointEdges.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingPointFaces.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingPointPoints.C
  ${polyMeshGenAddressing}/polyMeshGenAddressingParallelAddressing.C
  ${polyMeshGenChecks}/polyMeshGenChecks.C
  ${polyMeshGenChecks}/polyMeshGenChecksGeometry.C
  ${polyMeshGenChecks}/polyMeshGenChecksTopology.C
  ${partTetMesh}/partTetMesh.C
  ${partTetMesh}/partTetMeshAddressing.C
  ${partTetMesh}/partTetMeshParallelAddressing.C
  ${partTetMesh}/partTetMeshSimplex.C
  ${partTriMesh}/partTriMesh.C
  ${partTriMesh}/partTriMeshAddressing.C
  ${partTriMesh}/partTriMeshParallelAddressing.C
  ${partTriMesh}/partTriMeshSimplex.C
  ${triSurf}/triSurf.C
  ${triSurf}/triSurfPoints.C
  ${triSurf}/triSurfFacets.C
  ${triSurf}/triSurfFeatureEdges.C
  ${triSurf}/triSurfAddressing.C
  ${triSurf}/triSurfModifier.C
  ${findCellsIntersectingSurface}/findCellsIntersectingSurface.C
  ${tetCreatorOctree}/tetCreatorOctree.C
  ${tetCreatorOctree}/tetCreatorOctreePointsAndAddressing.C
  ${tetCreatorOctree}/tetCreatorOctreeFromFacesWithCentreNode.C
  ${tetCreatorOctree}/tetCreatorOctreeTetsAroundEdges.C
  ${tetCreatorOctree}/tetCreatorOctreeTetsAroundSplitEdges.C
  ${tetCreatorOctree}/tetCreatorOctreeTetsFromSplitFaces.C
  ${faceDecomposition}/faceDecomposition.C
  ${faceDecomposition}/decomposeFaces.C
  ${helperFunctions}/helperFunctionsStringConversion.C
  ${sortEdgesIntoChains}/sortEdgesIntoChains.C
  ${surfaceMorpherCells}/surfaceMorpherCells.C
  ${surfaceMorpherCells}/surfaceMorpherCellsMorphInternalFaces.C
  ${surfaceMorpherCells}/surfaceMorpherCellsCreateBoundaryFaces.C
  ${decomposeCells}/decomposeCells.C
  ${decomposeCells}/decomposeCellsPyramids.C
  ${decomposeCells}/decomposeCellsDecomposition.C
  ${topologicalCleaner}/topologicalCleaner.C
  ${topologicalCleaner}/topologyCleanerNonConsecutiveBoundaryVertices.C
  ${topologicalCleaner}/topologicalCleanerInvalidVertices.C
  ${topologicalCleaner}/topologyCleanerNonMappableCells.C
  ${checkCellConnectionsOverFaces}/checkCellConnectionsOverFaces.C
  ${checkIrregularSurfaceConnections}/checkIrregularSurfaceConnections.C
  ${checkIrregularSurfaceConnections}/checkIrregularSurfaceConnectionsFunctions.C
  ${checkNonMappableCellConnections}/checkNonMappableCellConnections.C
  ${checkBoundaryFacesSharingTwoEdges}/checkBoundaryFacesSharingTwoEdges.C
  ${boundaryLayers}/boundaryLayers.C
  ${boundaryLayers}/boundaryLayersCreateVertices.C
  ${boundaryLayers}/boundaryLayersFacesAndCells.C
  ${boundaryLayers}/boundaryLayerCells.C
  ${boundaryLayers}/boundaryLayersCheckTopologyOfBndFaces.C
  ${boundaryLayers}/boundaryLayersWrapperLayer.C
  ${extrudeLayer}/extrudeLayer.C
  ${detectBoundaryLayers}/detectBoundaryLayers.C
  ${detectBoundaryLayers}/detectBoundaryLayersFunctions.C
  ${refineBoundaryLayers}/refineBoundaryLayers.C
  ${refineBoundaryLayers}/refineBoundaryLayersFunctions.C
  ${refineBoundaryLayers}/refineBoundaryLayersFaces.C
  ${refineBoundaryLayers}/refineBoundaryLayersCells.C
  ${triangulateNonPlanarBaseFaces}/triangulateNonPlanarBaseFaces.C
  ${triangulateNonPlanarBaseFaces}/triangulateNonPlanarBaseFacesFunctions.C
  ${meshSurfaceEngine}/meshSurfaceEngine.C
  ${meshSurfaceEngine}/meshSurfaceEngineCalculateBoundaryNodesAndFaces.C
  ${meshSurfaceEngine}/meshSurfaceEngineParallelAddressing.C
  ${meshSurfaceEngine}/meshSurfaceEngineModifier.C
  ${meshSurfacePartitioner}/meshSurfacePartitioner.C
  ${meshSurfacePartitioner}/meshSurfacePartitionerFunctions.C
  ${meshSurfaceCheckInvertedVertices}/meshSurfaceCheckInvertedVertices.C
  ${meshSurfaceCheckEdgeTypes}/meshSurfaceCheckEdgeTypes.C
  ${meshSurfaceMapper}/meshSurfaceMapper.C
  ${meshSurfaceMapper}/meshSurfaceMapperMapVertices.C
  ${meshSurfaceMapper}/meshSurfaceMapperCornersAndEdges.C
  ${meshSurfaceMapper}/meshSurfaceMapperPremapVertices.C
  ${meshSurfaceMapper2D}/meshSurfaceMapper2D.C
  ${meshSurfaceMapper2D}/meshSurfaceMapper2DMapVertices.C
  ${meshSurfaceMapper2D}/meshSurfaceMapper2DPremapVertices.C
  ${edgeExtractor}/edgeExtractor.C
  ${edgeExtractor}/edgeExtractorCorners.C
  ${meshSurfaceEdgeExtractorNonTopo}/meshSurfaceEdgeExtractorNonTopo.C
  ${meshSurfaceEdgeExtractorNonTopo}/meshSurfaceEdgeExtractorNonTopoDistributeFaces.C
  ${meshSurfaceEdgeExtractor2D}/meshSurfaceEdgeExtractor2D.C
  ${meshSurfaceEdgeExtractor2D}/meshSurfaceEdgeExtractor2DDistributeFaces.C
  ${meshSurfaceEdgeExtractorFUN}/meshSurfaceEdgeExtractorFUN.C
  ${meshSurfaceEdgeExtractorFUN}/meshSurfaceEdgeExtractorFUNDistributeFaces.C
  ${createFundamentalSheets}/createFundamentalSheets.C
  ${createFundamentalSheetsFJ}/createFundamentalSheetsFJ.C
  ${createFundamentalSheetsJFS}/createFundamentalSheetsJFS.C
  ${correctEdgesBetweenPatches}/correctEdgesBetweenPatches.C
  ${correctEdgesBetweenPatches}/correctEdgesBetweenPatchesDistributeFaces.C
  ${renameBoundaryPatches}/renameBoundaryPatches.C
  ${meshOptimizer}/meshOptimizer.C
  ${meshOptimizer}/meshOptimizerOptimizePoint.C
  ${meshOptimizer}/meshOptimizerOptimizePointParallel.C
  ${meshOptimizer}/meshOptimizerOptimizeSurface.C
  ${meshOptimizer}/optimizeMeshFV.C
  ${tetMeshOptimisation}/tetMeshOptimisation.C
  ${tetMeshOptimisation}/tetMeshOptimisationParallel.C
  ${symmetryPlaneOptimisation}/symmetryPlaneOptimisation.C
  ${boundaryLayerOptimisation}/boundaryLayerOptimisation.C
  ${boundaryLayerOptimisation}/boundaryLayerOptimisationFunctions.C
  ${boundaryLayerOptimisation}/boundaryLayerOptimisationNormals.C
  ${boundaryLayerOptimisation}/boundaryLayerOptimisationThickness.C
  ${simplexSmoother}/simplexSmoother.C
  ${knuppMetric}/knuppMetric.C
  ${meshUntangler}/meshUntangler.C
  ${meshUntangler}/meshUntanglerCutRegion.C
  ${meshUntangler}/meshUntanglerCutRegionPoints.C
  ${meshUntangler}/meshUntanglerCutRegionEdges.C
  ${meshUntangler}/meshUntanglerCutRegionFaces.C
  ${meshUntangler}/meshUntanglerCutRegionTieBreak.C
  ${quadricMetric}/quadricMetric.C
  ${volumeOptimizer}/volumeOptimizer.C
  ${volumeOptimizer}/volumeOptimizerEvaluateGradients.C
  ${meshSurfaceOptimizer}/meshSurfaceOptimizer.C
  ${meshSurfaceOptimizer}/meshSurfaceOptimizerCalculateTrianglesAndAddressing.C
  ${meshSurfaceOptimizer}/meshSurfaceOptimizerOptimizePoint.C
  ${meshSurfaceOptimizer}/meshSurfaceOptimizerOptimizeSurface.C
  ${meshSurfaceOptimizer}/meshSurfaceOptimizerOptimizePointParallel.C
  ${surfaceOptimizer}/surfaceOptimizer.C
  ${meshOctreeCube}/meshOctreeCube.C
  ${meshOctreeCube}/meshOctreeCubeIntersections.C
  ${meshOctreeCube}/meshOctreeCubeRecursiveFunctions.C
  ${meshOctreeCube}/meshOctreeCubeRefine.C
  ${meshOctreeCube}/meshOctreeCubeCoordinatesIntersections.C
  ${meshOctreeModifier}/meshOctreeModifier.C
  ${meshOctreeModifier}/meshOctreeModifierRefineSelectedBoxes.C
  ${meshOctreeModifier}/meshOctreeModifierEnsureCorrectRegularity.C
  ${meshOctreeModifier}/meshOctreeModifierParallelRefinement.C
  ${meshOctreeModifier}/meshOctreeModifierDistributeLeavesToProcessors.C
  ${meshOctreeModifier}/meshOctreeModifierLoadDistribution.C
  ${meshOctreeModifier}/meshOctreeModifierReduceMemoryConsumption.C
  ${meshOctreeModifier}/meshOctreeModifierUpdateCommunicationPattern.C
  ${meshOctreeInsideOutside}/meshOctreeInsideOutside.C
  ${meshOctreeCreator}/meshOctreeCreator.C
  ${meshOctreeCreator}/meshOctreeCreatorAdjustOctreeToSurface.C
  ${meshOctreeCreator}/meshOctreeCreatorCreateOctreeBoxes.C
  ${meshOctreeCreator}/meshOctreeCreatorFrontalMarking.C
  ${meshOctreeCreator}/meshOctreeCreatorLoadDistribution.C
  ${meshOctreeAddressing}/meshOctreeAddressing.C
  ${meshOctreeAddressing}/meshOctreeAddressingCreation.C
  ${meshOctreeAddressing}/meshOctreeAddressingGluedMesh.C
  ${meshOctreeAddressing}/meshOctreeAddressingIrregularConnections.C
  ${meshOctreeAddressing}/meshOctreeAddressingParallelAddressing.C
  ${meshOctreeAutoRef}/meshOctreeAutomaticRefinement.C
  ${meshOctreeAutoRef}/meshOctreeAutomaticRefinementRef.C
  ${meshOctree}/meshOctree.C
  ${meshOctree}/meshOctreeCubePatches.C
  ${meshOctree}/meshOctreeNeighbourSearches.C
  ${meshOctree}/meshOctreeFindNearestSurfacePoint.C
  ${meshOctree}/meshOctreeInsideCalculations.C
  ${meshOctree}/meshOctreeParallelCommunication.C
  ${patchRefinement}/patchRefinement.C
  ${objectRefinement}/objectRefinement.C
  ${objectRefinement}/newObjectRefinement.C
  ${objectRefinement}/sphereRefinement.C
  ${objectRefinement}/lineRefinement.C
  ${objectRefinement}/coneRefinement.C
  ${objectRefinement}/boxRefinement.C
  ${objectRefinement}/hollowConeRefinement.C
  ${triSurface2DCheck}/triSurface2DCheck.C
  ${triSurfaceChecks}/triSurfaceChecks.C
  ${triSurfaceCleanupDuplicates}/triSurfaceCleanupDuplicates.C
  ${triSurfaceCleanupDuplicates}/triSurfaceCleanupDuplicatesFunctions.C
  ${triSurfaceCleanupDuplicateTriangles}/triSurfaceCleanupDuplicateTriangles.C
  ${triSurfaceCleanupDuplicateTriangles}/triSurfaceCleanupDuplicateTrianglesFunctions.C
  ${triSurfaceCopyParts}/triSurfaceCopyParts.C
  ${triSurfacePartitioner}/triSurfacePartitioner.C
  ${triSurfacePartitioner}/triSurfacePartitionerCreateAddressing.C
  ${triSurfaceCurvatureEstimator}/triSurfaceCurvatureEstimator.C
  ${triSurfaceCurvatureEstimator}/triSurfaceCurvatureEstimatorCalculate.C
  ${triSurfaceDetectFeatureEdges}/triSurfaceDetectFeatureEdges.C
  ${triSurfaceDetectFeatureEdges}/triSurfaceDetectFeatureEdgesFunctions.C
  ${triSurfaceClassifyEdges}/triSurfaceClassifyEdges.C
  ${triSurfaceClassifyEdges}/triSurfaceClassifyEdgesFunctions.C
  ${triSurfaceImportSurfaceAsSubset}/triSurfaceImportSurfaceAsSubset.C
  ${triSurfacePatchManipulator}/triSurfacePatchManipulator.C
  ${triSurfacePatchManipulator}/triSurfacePatchManipulatorFunctions.C
  ${triSurfaceRemoveFacets}/triSurfaceRemoveFacets.C
  ${triSurfaceRemoveFacets}/triSurfaceRemoveFacetsFunctions.C
  ${triSurfaceExtrude2DEdges}/triSurfaceExtrude2DEdges.C
  ${triSurfaceMetaData}/triSurfaceMetaData.C
  ${cartesianMeshExtractor}/cartesianMeshExtractor.C
  ${cartesianMeshExtractor}/cartesianMeshExtractorPointsAndAddressing.C
  ${cartesianMeshExtractor}/cartesianMeshExtractorPolyMesh.C
  ${cartesianMeshExtractor}/cartesianMeshExtractorDecomposeSplitHexes.C
  ${cartesianMeshGenerator}/cartesianMeshGenerator.C
  ${cartesian2DMeshGenerator}/cartesian2DMeshGenerator.C
  ${tetMeshExtractorOctree}/tetMeshExtractorOctree.C
  ${tetMeshGenerator}/tetMeshGenerator.C
  ${voronoiMeshExtractor}/voronoiMeshExtractor.C
  ${voronoiMeshExtractor}/voronoiMeshExtractorAddressing.C
  ${voronoiMeshExtractor}/voronoiMeshExtractorCreateMesh.C
  ${voronoiMeshGenerator}/voronoiMeshGenerator.C
  ${writeAsFPMA}/writeMeshFPMA.C
  ${writeAsFPMA}/fpmaMesh.C
  ${workflowControls}/workflowControls.C
)

add_foam_library(cfMesh SHARED ${SOURCES})
target_link_libraries(cfMesh PUBLIC meshTools edgeMesh)

if(OPENMP_FOUND)
  target_compile_definitions(cfMesh PUBLIC USE_OMP)
  target_compile_options(cfMesh PUBLIC ${OpenMP_CXX_FLAGS})
  target_link_libraries(cfMesh PUBLIC ${OpenMP_CXX_FLAGS})
endif()
