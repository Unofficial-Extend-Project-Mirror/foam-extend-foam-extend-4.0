//void Foam::edgeCorrectedVolPointInterpolation::makeExtrapolationWeights() const
//{
// Interpolation:
// For each point to be corrected calculate the extrapolated value
// for each face around it and combine them using the inverse
// distance weighting factors

//const labelList& ptc = boundaryPoints();

// Calculate the correction vectors

//extrapolationVectorsPtr_ =
//  new FieldField<Field, vector>(ptc.size());
//FieldField<Field, vector>& extraVecs = *extrapolationVectorsPtr_;
FieldField<Field, vector> extraVecs(ptc.size());

{
    const labelListList& pfaces = mesh.pointFaces();

    const volVectorField& centres = mesh.C();

    const fvBoundaryMesh& bm = mesh.boundary();

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = pfaces[curPoint];

        // extraVecs.hook(new vectorField(curFaces.size())); //- no hook function
        extraVecs.set
        (
            pointI,
            new vectorField(curFaces.size())
        );

        vectorField& curExtraVectors = extraVecs[pointI];

        label nFacesAroundPoint = 0;

        const vector& pointLoc = mesh.points()[curPoint];

        // Go through all the faces
        forAll (curFaces, faceI)
        {
            if (!mesh.isInternalFace(curFaces[faceI]))
            {
                // This is a boundary face.  If not in the empty patch
                // or coupled calculate the extrapolation vector
                label patchID =
                    mesh.boundaryMesh().whichPatch(curFaces[faceI]);

                if
                (
                    !isA<emptyFvPatch>(bm[patchID]) && !bm[patchID].coupled()
                )
                {
                    // Found a face for extrapolation
                    curExtraVectors[nFacesAroundPoint] =
                        pointLoc
                      - centres.boundaryField()[patchID]
                        [bm[patchID].patch().whichFace(curFaces[faceI])];

                    nFacesAroundPoint++;
                }
            }
        }

        curExtraVectors.setSize(nFacesAroundPoint);
    }
}
