//void volPointInterpolation::makeBoundaryWeights() const
//{
// const labelList& ptc = boundaryPoints();

// Calculate the correction vectors and weighting factors
//pointBoundaryWeightsPtr_ = new FieldField<Field, scalar>(ptc.size());
//FieldField<Field, scalar>& w = *pointBoundaryWeightsPtr_;
FieldField<Field, scalar> w(ptc.size());

{
    const labelListList& pf = mesh.pointFaces();

    const volVectorField& centres = mesh.C();

    const fvBoundaryMesh& bm = mesh.boundary();

    pointScalarField volPointSumWeights
    (
        IOobject
        (
            "volPointSumWeights",
            mesh.polyMesh::instance(),
            mesh
        ),
        pMesh,
        dimensionedScalar("zero", dimless, 0)
    );

    forAll (ptc, pointI)
    {
        const label curPoint = ptc[pointI];

        const labelList& curFaces = pf[curPoint];

        //w.hook(new scalarField(curFaces.size())); //philipc no hook function
        w.set
       (
           pointI,
           new scalarField(curFaces.size())
       );

       scalarField& curWeights = w[pointI];

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
                   !isA<emptyFvPatch>(bm[patchID])
                   && !(
                       bm[patchID].coupled()
                       //&& Pstream::parRun()
                       //&& !mesh.parallelData().cyclicParallel()
                   )
               )
               {
                   curWeights[nFacesAroundPoint] =
                       1.0/mag
                       (
                           pointLoc
                         - centres.boundaryField()[patchID]
                           [
                               bm[patchID].patch().whichFace(curFaces[faceI])
                           ]
                       );

                   nFacesAroundPoint++;
               }
           }
       }

       // Reset the sizes of the local weights
       curWeights.setSize(nFacesAroundPoint);

       // Collect the sum of weights for parallel correction
       volPointSumWeights[curPoint] += sum(curWeights);
    }

    // Do parallel correction of weights

    // Update coupled boundaries
    // Work-around for cyclic parallels.
    /*if (Pstream::parRun() && !mesh.parallelData().cyclicParallel())
    {
        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].initAddField();
            }
        }

        forAll (volPointSumWeights.boundaryField(), patchI)
        {
            if (volPointSumWeights.boundaryField()[patchI].coupled())
            {
                volPointSumWeights.boundaryField()[patchI].addField
                (
                    volPointSumWeights.internalField()
                );
            }
        }
    }*/

    // Re-scale the weights for the current point
    forAll (ptc, pointI)
    {
        w[pointI] /= volPointSumWeights[ptc[pointI]];
    }
}
