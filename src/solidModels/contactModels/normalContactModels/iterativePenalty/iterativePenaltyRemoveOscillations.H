   // oscillations sometimes appear in normal contact pressure
   // so I will try to limit them here
   // we will weight the current slavePressure_ the average of the
   // neighbours using the weight oscillationCorrectionFactor_
     //Pout << "Applying contact oscillation correction..." << flush;
{
  // slavePressure_

  const labelListList& faceFaces = slaveFaceZonePatch.faceFaces();

  // create global slavePressure
  vectorField globalSlavePressure(slaveFaceZonePatch.size(), vector::zero);
  const label slavePatchStart
    = mesh.boundaryMesh()[slavePatchIndex].start();
  forAll(slavePressure_, i)
    {
      globalSlavePressure
          [
              mesh.faceZones()[slaveFaceZoneID()
                  ].whichFace(slavePatchStart + i)] =
    slavePressure_[i];
    }
  // sum because each face is only on one proc
  reduce(globalSlavePressure, sumOp<vectorField>());

  // smooth slavePressure_ with face face disps
  forAll(faceFaces, facei)
    {
      label numFaceFaces = faceFaces[facei].size();

      if (mag(globalSlavePressure[facei]) > SMALL)
    //&&
    //numFaceFaces > 1) // don't smooth end/corner faces
    {
      vector avPress = vector::zero;
      int numNei = 0;
      forAll(faceFaces[facei], ffi)
        {
          label faceFace = faceFaces[facei][ffi];

          // only use faces in contact for average
          //if (slaveDispMag[faceFace] < -SMALL)
          //if (globalSlaveDispMag[faceFace] < -SMALL)
          // try averaging using all faces
          {
        avPress += globalSlavePressure[faceFace];
        numNei++;
          }
        }

      if (numNei == 0)
        {
          avPress = globalSlavePressure[facei];
        }
      else
        {
          avPress /= faceFaces[facei].size();
        }

         if (numFaceFaces == 1)
           {
         // for corner/end faces, decrease the weight of the neighbours
         avPress += globalSlavePressure[facei];
         avPress /= 2;
           }

         // weighted-average with face-faces
         globalSlavePressure[facei] =
           oscillationCorrFac_*globalSlavePressure[facei]
             + (1.0-oscillationCorrFac_)*avPress;
       }
       }

     // convert global back to local
     forAll(slavePressure_, facei)
       {
     slavePressure_[facei] =
       globalSlavePressure
       [
        mesh.faceZones()[slaveFaceZoneID()].whichFace(slavePatchStart + facei)
        ];
       }

     //Pout << "\tdone" << endl;
   }
