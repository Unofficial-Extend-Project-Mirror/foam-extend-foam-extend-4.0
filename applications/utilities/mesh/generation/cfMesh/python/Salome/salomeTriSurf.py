#!python
# =============================================================================
# Python module for writing OpenFOAM feature edge and triSurface files from
# within Salome platform.
# Tested on Salome 7.4.0 and python 2.7 on 64-bit Windows
#
# Author: Ivor Clifford <ivor.clifford@psi.ch>
#
# =============================================================================

def foamHeader(className, objectName):
    '''
    Return the OpenFOAM file header block as a string.
    '''
    return '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  2.2.1                                 |
|   \\\\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       %s;
    object      %s;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

''' % (className, objectName)


class triSurf:
    def __init__(self, object = None, allEdges = False):
        '''
        Construct from the supplied Salome mesh.

            object - the mesh object (must be a triangular surface mesh). If no
                    object is supplied, the current Salome selection is used.
            allEdges - If true, all edges on the mesh are included as feature
                    edges, otherwise only edges contained in groups are included

        NOTE: All face groups are assumed to represent patches. No face subsets
        are written. All edge groups are added as feature edge subsets. All point
        groups are added as point subsets.
        '''
        # =============================================================================
        # Initialize salome
        import salome
        import SMESH, SALOMEDS
        from salome.smesh import smeshBuilder
        from operator import itemgetter
        from collections import OrderedDict
        
        import SMESH, SALOMEDS
        
        # =============================================================================
        # Get the Salome mesh object
        if object is None:
            selected = salome.sg.getAllSelected()
            if len(selected) != 1:
                raise RuntimeError('A single Salome mesh object must be selected.')
            
            object = salome.myStudy.FindObjectID(selected[0]).GetObject()
        
        try:
            object.GetMesh()
        except:
            raise RuntimeError('Supplied object is not a Salome SMESH mesh.')
        
        smesh = smeshBuilder.New(salome.myStudy)
        mesh = smesh.Mesh(object)
        
        print "Converting SMESH Mesh '%s'" % mesh.GetName()
        
        # =============================================================================
        # Get basic mesh info
        nNodes = mesh.NbNodes()
        nFaces = mesh.NbFaces()
        nTris = mesh.NbTriangles()
        nEdges = mesh.NbEdges()
        nodeIds = mesh.GetNodesId()
        faceIds = mesh.GetElementsByType(SMESH.FACE)
        edgeIds = mesh.GetElementsByType(SMESH.EDGE)
        
        # Check that mesh is strictly triangular
        if nFaces != nTris:
            raise RuntimeError('Mesh is not strictly triangular')
        
        # Get patch and subset names & ids
        # All SMESH.FACE groups are assumed to be patches
        # All SMESH.EDGE groups are assumed to be feature subsets
        # All SMESH.NODE groups are assumed to be point subsets
        patches = OrderedDict()
        pointSubsets = OrderedDict()
        featureEdgeSubsets = OrderedDict()
        
        for group in mesh.GetGroups():
            if group.GetType() == SMESH.FACE:
                patches[group.GetName()] = group.GetIDs()
            elif group.GetType() == SMESH.EDGE:
                featureEdgeSubsets[group.GetName()] = group.GetIDs()
            elif group.GetType() == SMESH.NODE:
                pointSubsets[group.GetName()] = group.GetIDs()
        
        # =============================================================================
        # Process faces and patches
        # Get patchId for each face
        lastPatchId = len(patches)
        patchIds = [lastPatchId] * max(faceIds)
        patchId = 0
        for name, ids in patches.iteritems():
            for faceId in ids:
                if patchIds[faceId-1] == lastPatchId:
                    patchIds[faceId-1] = patchId
                else:
                    print "Face %d is assigned to both groups %s and %s" % (faceId, name, patches.keys()[patchIds[faceId-1]])
                    raise RuntimeError('Groups of faces are not unique, i.e. they overlap.')
                
            patchId += 1
        
        # Compact and reorder patchIds to match faceIds
        patchIds = [patchIds[faceId-1] for faceId in faceIds]
        
        # Reorder faces by increasing group id
        faceAndpatchIds = sorted(zip(faceIds, patchIds), key=itemgetter(1))
        faceIds, patchIds = zip(*faceAndpatchIds)
        
        # Add unused faces to the default patch
        defaultFaces = [faceId for faceId, patchId in faceAndpatchIds if patchId == lastPatchId]
        if len(defaultFaces) > 0:
            patches['defaultFaces'] = defaultFaces
        
        defaultFaces = None
        
        # =============================================================================
        # Process feature edges
        if not allEdges:
            edgeIds = []
            for name, ids in featureEdgeSubsets.iteritems():
                edgeIds += ids
            
            edgeIds = list(set(edgeIds))
            nEdges = len(edgeIds)
        
        # Reverse mapping of edge ids since they aren't necessarily numbered 1..nEdges
        if len(edgeIds):
            edgeMap = [-1] * max(edgeIds)
        else:
            edgeMap = []
        
        i=0
        for edgeId in edgeIds:
            edgeMap[edgeId-1] = i
            i += 1
        
        # =============================================================================
        # Process nodes
        # Reverse mapping of node ids since nodes aren't necessarily numbered 1..nNodes
        nodeMap = [-1] * max(nodeIds)
        i=0
        for nodeId in nodeIds:
            nodeMap[nodeId-1] = i
            i += 1
        
        # =============================================================================

        self._mesh = mesh
        
        self._nodeIds = nodeIds
        self._edgeIds = edgeIds
        self._faceIds = faceIds
        
        self._nodeMap = nodeMap
        self._edgeMap = edgeMap
        self._faceMap = []
        
        self._patches = patches
        self._pointSubsets = pointSubsets
        self._featureEdgeSubsets = featureEdgeSubsets
        self._faceSubsets = {}
        
        print 'Done'
    
    def nNodes(self):
        '''
        Return the number of nodes
        '''
        return len(self._nodeIds)
    
    def nEdges(self):
        '''
        Return the number of edges
        '''
        return len(self._edgeIds)
    
    def nFacets(self):
        '''
        Return the number of triangular facets
        '''
        return len(self._faceIds)
    
    def nPatches(self):
        '''
        Return the number of patches
        '''
        return len(self._patches)
    
    def _writePatchDefs(self, f, typeName = 'wall'):
        '''
        Write the patch definitions to file as an OpenFOAM geometricSurfacePatchList.
        NOTE: All patches are assumed to be walls.
        '''
        patches = self._patches
        
        f.write('%d\n(\n' % len(patches))
        for name in patches.iterkeys():
            f.write('%s\t%s\n' % (name, typeName))
        
        f.write(')\n')
    
    def _writeNodes(self, f):
        '''
        Write the nodes to file as an OpenFOAM pointField.
        '''
        mesh = self._mesh
        nodeIds = self._nodeIds
        
        f.write('%d\n(\n' % len(nodeIds))
        
        for x, y, z in [mesh.GetNodeXYZ(nodeId) for nodeId in nodeIds]:
            f.write( '( %g %g %g )\n' % (x, y, z))
        
        f.write(')\n')
    
    def _writeFeatureEdges(self, f):
        '''
        Write the feature edges to file as an OpenFOAM edgeList.
        '''
        mesh = self._mesh
        nodeMap = self._nodeMap
        edgeIds = self._edgeIds
        
        f.write('%d\n(\n' % len(edgeIds))
        
        for edgeId in edgeIds:
            nodes = mesh.GetElemNodes(edgeId)
            f.write( '(' + ' '.join([str(nodeMap[nodeId-1]) for nodeId in nodes]) + ')\n')
        
        f.write(')\n')
    
    def _writeFacets(self, f):
        '''
        Write the facets to file as an OpenFOAM List of labelledTri.
        '''
        from itertools import chain
        
        mesh = self._mesh
        nodeMap = self._nodeMap
        patches = self._patches
        
        f.write('%d\n(\n' % sum([len(patch) for patch in patches.itervalues()]))
        
        patchId = 0
        for patchId, (patchName, faceIds) in enumerate(patches.iteritems()):
            for faceId in faceIds:
                nodes = mesh.GetElemNodes(faceId)
                f.write( '((' + ' '.join([str(nodeMap[nodeId-1]) for nodeId in nodes]) + ') %d)\n' % patchId)
        
        f.write(')\n')
    
    def _writeSubsets(self, f, subsets, map, typeId):
        '''
        General function to write a subset to file as an OpenFOAM Map<meshSubset>.
        '''
        f.write('%d\n(\n' % len(subsets))
        for name, ids in subsets.iteritems():
            f.write('%s %s %d ( %s )\n' % (name, typeId, len(ids), ' '.join([str(map[id-1]) for id in ids])))
        
        f.write(')\n')
    
    def _writePointSubsets(self, f):
        '''
        Write the point subsets to file as and OpenFOAM Map<meshSubset>.
        '''
        self._writeSubsets(f, self._pointSubsets, self._nodeMap, '2')
    
    def _writeFaceSubsets(self, f):
        '''
        Write the face subsets to file as and OpenFOAM Map<meshSubset>.
        '''
        self._writeSubsets(f, self._faceSubsets, self._faceMap, '4')
        
    def _writeFeatureEdgeSubsets(self, f):
        '''
        Write the feature edge subsets to file as and OpenFOAM Map<meshSubset>.
        '''
        self._writeSubsets(f, self._featureEdgeSubsets, self._edgeMap, '8')
    
    def writeEdgeMesh(self, fileName):
        '''
        Write to file as an OpenFOAM edgeMesh
        
            fileName - The file name to write
        '''
        # Create file
        f = open(fileName, 'wb') # NOTE: file opened as binary to ensure unix-style line breaks
        
        # Write header
        f.write(foamHeader('edgeMesh', self._mesh.GetName()))
        
        self._writeNodes(f)
        self._writeFeatureEdges(f)
        
        f.close()
        
        print 'edgeMesh written to %s' % fileName
    
    def writeFtr(self, fileName):
        '''
        Write to file as an OpenFOAM cfMesh FTR file
        
            fileName - the file name to write
        '''
        # Create file
        f = open(fileName, 'wb') # NOTE: file opened as binary to ensure unix-style line breaks
        
        self._writePatchDefs(f)
        self._writeNodes(f)
        self._writeFacets(f)
        
        f.close()
        
        print 'triSurf written to %s' % fileName
        
    def writeFms(self, fileName):
        '''
        Write to file as an OpenFOAM cfMesh FMS file
        
            fileName - the file name to write
        '''
        # Create file
        f = open(fileName, 'wb') # NOTE: file opened as binary to ensure unix-style line breaks
        
        self._writePatchDefs(f)
        self._writeNodes(f)
        self._writeFacets(f)
        self._writeFeatureEdges(f)
        self._writePointSubsets(f)
        self._writeFaceSubsets(f)
        self._writeFeatureEdgeSubsets(f)
        
        f.close()
        
        print 'triSurf written to %s' % fileName


if __name__ == '__main__':
    import salome
    salome.salome_init()
    
    import SMESH, SALOMEDS
