#!python
# =============================================================================
# Salome GEOM script to extract the feature edges from a body and add them
# to the group named 'featureEdges'.
# Tested on Salome 7.4.0 and python 2.7 on 64-bit Windows
#
# Author: Ivor Clifford <ivor.clifford@psi.ch>
#
# =============================================================================

def extractFeatureEdges(body, minFeatureAngle = 5):
    '''
    Find all feature edges on the supplied body and return them as a list
    of edge ids.

        body - A Salome solid, compound, shell or face object to find all
                feature edges on.
        minFeatureAngle - the angle (in degrees) between adjacent surfaces
                above which the edge will be considered a feature angle.
    '''
    import salome
    salome.salome_init()

    import GEOM
    from salome.geom import geomBuilder
    geompy = geomBuilder.New(salome.myStudy)
    
    # Check the body type
    if not (body.GetShapeType() in [GEOM.SHELL, GEOM.SOLID, GEOM.FACE, GEOM.COMPOUND]):
       raise RuntimeError('Supplied object is not a solid, shell or face.')
    
    print 'Extracting edges of %s with feature angle > %g.' % (body.GetName(), minFeatureAngle)
    
    # Extract basic info
    faces = geompy.SubShapeAll(body, geompy.ShapeType["FACE"])
    curves = geompy.SubShapeAll(body, geompy.ShapeType["EDGE"])
    points = geompy.SubShapeAll(body, geompy.ShapeType["VERTEX"])
    
    faceIds = geompy.GetSubShapesIDs(body, faces)
    curveIds = geompy.GetSubShapesIDs(body, curves)
    nodeIds = geompy.GetSubShapesIDs(body, points)
    
    maxFaceId = max(faceIds)
    maxCurveId = max(curveIds)
    maxNodeId = max(nodeIds)
    
    # Reverse mapping from curve id to local curve arrays
    faceMap = [-1 for i in xrange(maxFaceId+1)]
    for localId, id in enumerate(faceIds):
        faceMap[id] = localId
    
    curveMap = [-1 for i in xrange(maxCurveId+1)]
    for localId, id in enumerate(curveIds):
        curveMap[id] = localId
    
    nodeMap = [-1 for i in xrange(maxNodeId+1)]
    for localId, id in enumerate(nodeIds):
        nodeMap[id] = localId
    
    
    # Get curves on each face
    faceCurveIds = [[curveMap[id] for id in geompy.GetSubShapesIDs(
            body,
            geompy.SubShapeAll(face, geompy.ShapeType["EDGE"])
        )] for face in faces]
    
    # Get faces attached to each curve
    curveFaceIds = [[] for id in curveIds]
    
    for faceI, ids in enumerate(faceCurveIds):
        for id in ids:
            curveFaceIds[id].append(faceI)
    
    # Now that we have the connectivity for curves and faces find the
    # feature edges
    featureEdgeIds = []
    for curveId, curve, adjFaceIds in zip(curveIds, curves, curveFaceIds):
        if len(adjFaceIds) == 2:
            # Curve with 2 adjacent faces - Test feature angle
            # Determine break angle at each curve
            # If greater than the feature edge angle, add the curve to group featureEdges
            face1 = faces[adjFaceIds[0]]
            face2 = faces[adjFaceIds[1]]
            point = geompy.GetFirstVertex(curve)    # Test at the first vertex
            n1 = geompy.GetNormal(face1, point)
            n2 = geompy.GetNormal(face2, point)
            angle = geompy.GetAngle(n1, n2)
            if angle > minFeatureAngle:
                featureEdgeIds.append(curveId)
        
        elif len(adjFaceIds) == 1:
            # Curve on standalone face - Add by default
            featureEdgeIds.append(curveId)
        
        elif len(adjFaceIds) == 0:
            # Standalone curve - Ignore
            None
        
        else:
            raise RuntimeError('Curve found sharing %d faces. This is unexpected for fully enclosed bodies.' % len(adjFaceIds))
    
    # Done
    print "%d feature edges found" % len(featureEdgeIds)
    
    return featureEdgeIds


# If run as a standalone script, use the current Salome GUI selection
# and add the feature edges to group named 'featureEdges'
if __name__ == '__main__':
    import salome
    salome.salome_init()
    
    import GEOM
    from salome.geom import geomBuilder
    geompy = geomBuilder.New(salome.myStudy)
    
    # Get current GUI selection
    selected = salome.sg.getAllSelected()
    if len(selected) != 1:
        raise RuntimeError('A single solid, shell or face object must be selected.')
    
    body = salome.myStudy.FindObjectID(selected[0]).GetObject()
    
    # Get feature edges and add to the group 'featureEdges'
    featureEdges = geompy.CreateGroup(body, geompy.ShapeType["EDGE"])
    geompy.UnionIDs(featureEdges, extractFeatureEdges(body))
    geompy.addToStudyInFather(body, featureEdges, 'featureEdges')
    
    if salome.sg.hasDesktop():
        salome.sg.updateObjBrowser(1)
