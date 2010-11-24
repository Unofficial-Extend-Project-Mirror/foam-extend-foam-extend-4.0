// Boundary processor patch points
{
    const labelList& curPointEdges = pointEdges[curPoint];

    label patchID = -1;
    label edgeID = -1;

    forAll (curPointEdges, edgeI)
    {
        label curEdge = curPointEdges[edgeI];

        if (edgeFaces[curEdge].size() == 1)
        {
            forAll(aMesh().boundary(), patchI)
            {
                const labelList& curEdges = aMesh().boundary()[patchI];

                label index = findIndex(curEdges, curEdge);

                if (index != -1)
                {
                    if 
                    (
                        aMesh().boundary()[patchI].type()
                     != processorFaPatch::typeName
                    )
                    {
                        patchID = patchI;
                        edgeID = index;
                        break;
                    }
                }
            }
        }
    }

    if (patchID != -1)
    {
        vector mirrorCtrlPoint =
            patchMirrorPoints[patchID][edgeID];

        lsPoints[curPatchPoint].setSize(lsPoints[curPatchPoint].size() + 1);
        lsPoints[curPatchPoint][lsPoints[curPatchPoint].size() - 1] =
            mirrorCtrlPoint;
    }
}
