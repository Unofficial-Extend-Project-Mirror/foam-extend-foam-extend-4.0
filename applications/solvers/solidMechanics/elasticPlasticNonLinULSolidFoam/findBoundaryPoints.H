//void volPointInterpolation::makeBoundaryAddressing() const
//{

// Go through all patches and mark up the external edge points
labelHashSet pointsCorrectionMap(max(mesh.nPoints()/10, 100));

const fvBoundaryMesh& bm = mesh.boundary();

forAll (bm, patchI)
{
    // If the patch is empty, skip it
    // If the patch is coupled, and there are no cyclic parallels, skip it
    if
    (
        !isA<emptyFvPatch>(bm[patchI])
        && !(
            bm[patchI].coupled()
            //&& Pstream::parRun()
            //&& !mesh.parallelData().cyclicParallel()
        )
    )
    {
        const labelList& bp = bm[patchI].patch().boundaryPoints();

        const labelList& meshPoints = bm[patchI].patch().meshPoints();

        forAll (bp, pointI)
        {
            pointsCorrectionMap.insert(meshPoints[bp[pointI]]);
        }
    }
}

// Grab the points to correct
//boundaryPointsPtr_ = new labelList(pointsCorrectionMap.toc());
//labelList& ptc = *boundaryPointsPtr_;
labelList ptc(pointsCorrectionMap.toc());

sort(ptc);
