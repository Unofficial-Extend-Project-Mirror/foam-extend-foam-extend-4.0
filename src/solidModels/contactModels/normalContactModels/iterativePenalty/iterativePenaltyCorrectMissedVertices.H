// sometimes the pointDistanceToIntersection "misses" vertices which
// are in contact so the pointGap is GREAT
// these missed vertices can cause serious problems with the contact convergence
// Here I will check for these missed vertices and set their contact gap
// value to the average of their non-GREAT neighbours
// essentially this is a quick fix but it seems to work nicely
  {
    const labelListList& pointPoints = *slavePointPointsPtr_;
    scalarField newGap = slavePointPenetration_;
    scalar bbFrac = 0.5*mesh.bounds().mag();
    forAll(slavePointPenetration_, pointi)
      {
    // for points with a contactGap of GREAT, we will
    // reset their gap to the average of their non-GREAT
    // neighbours
    if (slavePointPenetration_[pointi] > 0.9*GREAT)
      {
        numCorrectedPoints++;
        scalar averageNei = 0.0;
        label numNei = 0;
        forAll(pointPoints[pointi], ppi)
          {
        if (slavePointPenetration_[pointPoints[pointi][ppi]] < 0.9*GREAT)
          {
            averageNei += slavePointPenetration_[pointPoints[pointi][ppi]];
            numNei++;
          }
          }
        if (numNei > SMALL)
          {
        averageNei /= numNei;
          }
        else
          {
        // limit max gap to fraction of mesh bounding box
        averageNei = bbFrac;
          }

        newGap[pointi] = averageNei;
      }
      }

    slavePointPenetration_ = newGap;
  }
